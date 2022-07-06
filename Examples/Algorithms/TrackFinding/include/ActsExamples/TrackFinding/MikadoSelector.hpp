// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <vector>
#include <limits>

namespace ActsExamples {

class MikadoSelector {
public:
  using Config = Acts::MeasurementSelector::Config;

  MikadoSelector() = default;
  MikadoSelector(Config config, std::size_t n) 
    : m_config(config),
      m_used_meas(n, false)
  {}

  bool 
  hasBeenUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const {
    const auto& sourceLink =
      static_cast<const IndexSourceLink&>(trackstate.uncalibrated());
    std::size_t idMeas = sourceLink.index();
    return m_used_meas[idMeas];
  }

  Acts::Result<std::pair<std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator,
			 std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator>>
  select(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
	 bool& isOutlier, Acts::LoggerWrapper logger) const
  {
    std::cout<<">>> Running Mikado selector"<<std::endl;

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

    double minChi2 = std::numeric_limits<double>::max();
    std::size_t n_used_meas = 0;
    size_t minIndex = 0;
    size_t index = 0;
    // Loop over all measurements to select the compatible measurements
    for (auto& trackState : candidates) {      
      // Take the parameter covariance
      const auto predicted = trackState.predicted();
      const auto predictedCovariance = trackState.predictedCovariance();
      Acts::visit_measurement(
			trackState.calibrated(), trackState.calibratedCovariance(),
			trackState.calibratedSize(),
			[&](const auto calibrated, const auto calibratedCovariance) {
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

	  if (hasBeenUsed(trackState)) {
	    n_used_meas++;
	  } else if (chi2 < minChi2) {
	    // Search for the measurement with the min chi2
	    minChi2 = chi2;
	    minIndex = index;
	  }
        }); // visit measurement
      
      index++;
    }

    const auto& chi2CutOff = cuts->chi2CutOff;

    {
      // If there is no selected measurement, return the measurement with the best
      // chi2 and tag it as an outlier
      const auto bestIt = std::next(candidates.begin(), minIndex);
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

    std::sort(
	      candidates.begin(), candidates.end(),
	      [this](const auto& tsa, const auto& tsb) 
	      {
		if (hasBeenUsed(tsa)) return false;
		if (hasBeenUsed(tsb)) return true;
		return tsa.chi2() < tsb.chi2(); 
	      });

    // checking sorting
    // std::cout<<"Total used Measurements: "<< n_used_meas << std::endl;
    // std::cout<<"Checking sorting for " << candidates.size() << " candidates: "<<std::endl;
    // for (const auto& el : candidates) {
    //   bool used = hasBeenUsed(el);
    //   const auto chi2 = el.chi2();
    //   std::cout <<"   -- chi2=" << chi2 << " ; used=" << (used?"Y":"N") << std::endl;
    // }
    // std::cout << "done"<<std::endl;

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

    // mark as used()
    for (auto it = candidates.begin(); it!= endIterator; it++) {
      const auto& sourceLink =
	static_cast<const IndexSourceLink&>(it->uncalibrated());
      std::size_t idMeas = sourceLink.index();
      if (m_used_meas[idMeas]) std::cout<<"MARKING ALREADY TRUE MEASUREMENT IN MIKADO SELECTION" << std::endl;
      m_used_meas[idMeas] = true;
    }
    return std::pair{candidates.begin(), endIterator};
  }
  
  
private:
  template <typename cut_value_t>
  static cut_value_t VariableCut(
    const Acts::MultiTrajectory::TrackStateProxy& trackState,
    const Acts::MeasurementSelector::Config::Iterator selector,
    const std::vector<cut_value_t>& cuts, Acts::LoggerWrapper logger);
  
private:
  Config m_config;
  mutable std::vector<bool> m_used_meas;
};

template <typename cut_value_t>
cut_value_t MikadoSelector::VariableCut(
  const Acts::MultiTrajectory::TrackStateProxy& trackState,
  const Acts::MeasurementSelector::Config::Iterator selector,
  const std::vector<cut_value_t>& cuts, Acts::LoggerWrapper logger) 
{
  const auto& etaBins = selector->etaBins;
  if (etaBins.empty()) {
    return cuts[0];  // shortcut if no etaBins
  }
  const auto eta = std::atanh(std::cos(trackState.predicted()[Acts::eBoundTheta]));
  const auto abseta = std::abs(eta);
  size_t bin = 0;
  for (auto etaBin : etaBins) {
    if (etaBin >= abseta) {
      break;
    }
    bin++;
  }
  if (bin >= cuts.size()) {
    bin = cuts.size() - 1;
  }
  ACTS_VERBOSE("Variable cut for eta=" << eta << ": " << cuts[bin]);
  return cuts[bin];
}
  
} // namespace
