// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "ActsExamples/TrackFinding/MikadoSelector.hpp" 

namespace ActsExamples {

  Acts::Result<std::pair<std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator,
                         std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator>>
  MikadoSelector::select(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
			 bool& isOutlier, Acts::LoggerWrapper logger) const 
  {
    using Result = Acts::Result<
      std::pair<std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator,
		std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator>>;
    
    ACTS_VERBOSE("Invoked Mikado Selector");

    // Return error if no measurement
    if (candidates.empty()) {
      return Acts::CombinatorialKalmanFilterError::MeasurementSelectionFailed;
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
      return Acts::CombinatorialKalmanFilterError::MeasurementSelectionFailed;
    }

    // select best chi2
    std::size_t n_used_meas = 0;

    auto select_trackstate
      = [this, &n_used_meas] (const Acts::MultiTrajectory::TrackStateProxy& trackState) -> bool
      {
	bool has_been_used = hasBeenUsed(trackState);
	if (has_been_used) n_used_meas++;
	return not has_been_used;
      };
    const auto bestIt = selectBestCandidate(candidates, select_trackstate);

    
    const auto& chi2CutOff = cuts->chi2CutOff;

    {
      // If there is no selected measurement, return the measurement with the best
      // chi2 and tag it as an outlier
      const auto chi2 = bestIt->chi2();
      const auto chi2Cut = VariableCut(*bestIt, cuts, chi2CutOff, logger);
      ACTS_VERBOSE("Chi2: " << chi2 << ", max: " << chi2Cut);
      if (chi2 >= chi2Cut or n_used_meas >= candidates.size()) {
	ACTS_VERBOSE("No measurement candidate. Return an outlier measurement.");
	isOutlier = true;
	// return single item range, no sorting necessary
	return Result::success(std::pair{bestIt, std::next(bestIt, 1)});
      }
    }

    std::sort(candidates.begin(), candidates.end(),
	      [this] (const auto& tsa, const auto& tsb) -> bool
	      {
		if (hasBeenUsed(tsa)) return false;
		if (hasBeenUsed(tsb)) return true;
		return tsa.chi2() < tsb.chi2(); 
	      });


    // use |eta| of best measurement to select numMeasurementsCut
    const std::size_t numMeasurementsCut = std::min(
						    VariableCut(*candidates.begin(), cuts, cuts->numMeasurementsCutOff, logger),
						    candidates.size() - n_used_meas
						    );
    
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

    // Set as used
    for (auto it = candidates.begin(); it!= endIterator; it++) 
      setAsUsed(*it);

    ACTS_VERBOSE("Number of selected measurements: "
		 << std::distance(candidates.begin(), endIterator)
		 << ", max: " << numMeasurementsCut);

    isOutlier = false;

    return std::pair{candidates.begin(), endIterator};
  }
  
} // namespace
