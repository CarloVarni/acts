// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Acts {

Result<std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                 std::vector<MultiTrajectory::TrackStateProxy>::iterator>>
MeasurementSelector::select(
    std::vector<MultiTrajectory::TrackStateProxy>& candidates, bool& isOutlier,
    LoggerWrapper logger) const {
  using Result = Result<
      std::pair<std::vector<MultiTrajectory::TrackStateProxy>::iterator,
                std::vector<MultiTrajectory::TrackStateProxy>::iterator>>;

  ACTS_VERBOSE("Invoked MeasurementSelector");

  // Return error if no measurement
  if (candidates.empty()) {
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  // Get geoID of this surface
  auto surface = &candidates.front().referenceSurface();
  auto geoID = surface->geometryId();

  // Find the appropriate cuts
  auto cuts = m_config.find(geoID);
  if (cuts == m_config.end()) {
    // for now we consider missing cuts an unrecoverable error
    // TODO consider other options e.g. do not add measurements at all (not
    // even as outliers)
    return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
  }

  const auto bestIt = selectBestCandidate(candidates);

  const auto& chi2CutOff = cuts->chi2CutOff;

  {
    // If there is no selected measurement, return the measurement with the best
    // chi2 and tag it as an outlier
    const auto chi2 = bestIt->chi2();
    const auto chi2Cut = VariableCut(*bestIt, cuts, chi2CutOff, logger);
    ACTS_VERBOSE("Chi2: " << chi2 << ", max: " << chi2Cut);
    if (chi2 >= chi2Cut) {
      ACTS_VERBOSE("No measurement candidate. Return an outlier measurement.");
      isOutlier = true;
      // return single item range, no sorting necessary
      return Result::success(std::pair{bestIt, std::next(bestIt, 1)});
    }
  }

  std::sort(
      candidates.begin(), candidates.end(),
      [](const auto& tsa, const auto& tsb) { return tsa.chi2() < tsb.chi2(); });

  // use |eta| of best measurement to select numMeasurementsCut
  const auto numMeasurementsCut = VariableCut(
      *candidates.begin(), cuts, cuts->numMeasurementsCutOff, logger);

  auto endIterator = candidates.begin();
  auto maxIterator = candidates.end();
  if (candidates.size() > numMeasurementsCut && numMeasurementsCut > 0) {
    maxIterator = std::next(candidates.begin(), numMeasurementsCut);
  }

  ++endIterator;  // best measurement already confirmed good
  for (; endIterator != maxIterator; ++endIterator) {
    const auto chi2 = endIterator->chi2();
    const auto chi2Cut = VariableCut(*endIterator, cuts, chi2CutOff, logger);
    ACTS_VERBOSE("Chi2: " << chi2 << ", max: " << chi2Cut);
    if (chi2 >= chi2Cut) {
      break;  // endIterator now points at the first track state with chi2
              // larger than our cutoff => defines the end of our returned
              // range
    }
  }
  ACTS_VERBOSE("Number of selected measurements: "
               << std::distance(candidates.begin(), endIterator)
               << ", max: " << numMeasurementsCut);

  isOutlier = false;
  return std::pair{candidates.begin(), endIterator};
}

std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator	
MeasurementSelector::selectBestCandidate(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
					 std::function<bool(const Acts::MultiTrajectory::TrackStateProxy&)> selection_function) const
{
  double minChi2 = std::numeric_limits<double>::max();
  std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator it = candidates.begin();
  
  for (auto iter = candidates.begin(); iter != candidates.end(); iter++) {
    auto& trackState = *iter;
    const auto predicted = trackState.predicted();
    const auto predictedCovariance = trackState.predictedCovariance();
    
    if ( not selection_function(trackState) )
      continue;
    
    Acts::visit_measurement(
		       trackState.calibrated(), trackState.calibratedCovariance(),
		       trackState.calibratedSize(),
		       [&] (const auto calibrated, const auto calibratedCovariance) {
	      constexpr size_t kMeasurementSize =
		decltype(calibrated)::RowsAtCompileTime;
	      
	      using ParametersVector = Acts::ActsVector<kMeasurementSize>;
	      
	      // Take the projector (measurement mapping function)
	      const auto H =
		trackState.projector()
		.template topLeftCorner<kMeasurementSize, Acts::eBoundSize>()
		.eval();
	      
	      // Get the residuals
	      ParametersVector res;
	      res = calibrated - H * predicted;
	      
	      // Get the chi2
	      double& chi2 = trackState.chi2();
	      chi2 = (res.transpose() *
		      ((calibratedCovariance +
			H * predictedCovariance * H.transpose()))
		      .inverse() *
		      res)
		.eval()(0, 0);
	      
	      if (chi2 < minChi2) {
		// Search for the measurement with the min chi2
		minChi2 = chi2;
		it = iter;
	      }
	    }); // visit measurement
    
  } // for loop

  return it;
}
  
}  // namespace Acts
