// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

  template <typename traj_t>
  Result<std::pair<
	   typename std::vector<typename traj_t::TrackStateProxy>::iterator,
	   typename std::vector<typename traj_t::TrackStateProxy>::iterator>>
  MeasurementSelector::select(std::vector<typename traj_t::TrackStateProxy>& candidates,
			      bool& isOutlier, const Logger& logger) const
  {
    using Result = Result<std::pair<
      typename std::vector<typename traj_t::TrackStateProxy>::iterator,
			    typename std::vector<typename traj_t::TrackStateProxy>::iterator>>;
    auto start = std::chrono::high_resolution_clock::now();
    
    ACTS_VERBOSE("Invoked MeasurementSelector");

    // Return error if no measurement
    if (candidates.empty()) {
      return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
    }

    // Get geoID of this surface
    auto geoID = candidates.front().referenceSurface().geometryId();

    // Find the appropriate cuts
    auto cuts = m_config.find(geoID);
    if (cuts == m_config.end()) {
      // for now we consider missing cuts an unrecoverable error
      // TODO consider other options e.g. do not add measurements at all (not
      // even as outliers)
      return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
    }

    isOutlier = false;
    
    assert(!cuts->chi2CutOff.empty());
    const std::vector<double>& chi2CutOff = cuts->chi2CutOff;
    double maxChi2Cut = *std::max_element(chi2CutOff.begin(), chi2CutOff.end());

    double minChi2 = std::numeric_limits<double>::max();
    std::size_t minIndex = 0;

    // Loop over all measurements to select the compatible measurements
    // Sort track states which do not satisfy the chi2 cut to the end.
    // When done trackStateIterEnd will point to the first element that
    // does not satisfy the chi2 cut.    
    std::size_t goodMeasurements = 0ul;
    for (std::size_t i(0ul); i<candidates.size(); ++i) {
      auto& trackState = candidates[i];

      double chi2 = calculateChi2(
            // This abuses an incorrectly sized vector / matrix to access the
            // data pointer! This works (don't use the matrix as is!), but be
            // careful!
            trackState.template calibrated<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackState.template calibratedCovariance<
                    MultiTrajectoryTraits::MeasurementSizeMax>()
                .data(),
            trackState.predicted(), trackState.predictedCovariance(),
            trackState.projector(), trackState.calibratedSize());

      trackState.chi2() = chi2;

      if (chi2 < maxChi2Cut &&
	  chi2 < VariableCut<traj_t>(trackState, cuts, chi2CutOff, logger)) {
	if (goodMeasurements != i)
	  std::swap(candidates[goodMeasurements], candidates[i]);
	++goodMeasurements;
      }

      // get the best one since we need it anyway later
      if (chi2 < minChi2) {
	minChi2 = chi2;
	minIndex = i;
      }
    }
    
    const std::size_t numMeasurementsCut = VariableCut<traj_t>(candidates[minIndex],
							       cuts,
							       cuts->numMeasurementsCutOff,
							       logger);
    // If there are no measurements below the chi2 cut off, return the
    // measurement with the best chi2 and tag it as an outlier
    // The same applies if we have at least 1 measurement and the numMeasurementsCut is 1
    if (goodMeasurements == 0ul || numMeasurementsCut == 1ul) {
      const auto bestIt = candidates.begin() + minIndex;
      ACTS_VERBOSE(
          "No measurement candidate. Return an outlier measurement chi2="
          << bestIt->chi2());
      isOutlier = true;
      // return single item range, no sorting necessary
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
      std::cout << "Duration: " << duration << std::endl;
      return Result::success(std::make_pair(bestIt, bestIt+1));
    }

    if (goodMeasurements <= numMeasurementsCut) {
      ACTS_VERBOSE("Number of selected measurements: " << goodMeasurements
		   << ", max: " << numMeasurementsCut);
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
      std::cout << "Duration: " << duration << std::endl;
      return Result::success(std::make_pair(candidates.begin(), candidates.begin() + goodMeasurements));
    }
    
    std::sort(candidates.begin(), candidates.begin() + goodMeasurements,
              [] (const auto& tsa, const auto& tsb) {
                return tsa.chi2() < tsb.chi2();
              });

    ACTS_VERBOSE("Number of selected measurements: " << numMeasurementsCut
                 << ", max: " << numMeasurementsCut);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    std::cout << "Duration: " << duration << std::endl;
    return Result::success(std::make_pair(candidates.begin(), candidates.begin() + numMeasurementsCut));
  }

  template <typename traj_t, typename cut_value_t>
  cut_value_t MeasurementSelector::VariableCut(
      const typename traj_t::TrackStateProxy& trackState,
      const Acts::MeasurementSelector::Config::Iterator selector,
      const std::vector<cut_value_t>& cuts, const Logger& logger) {
    const auto& etaBins = selector->etaBins;
    if (etaBins.empty()) {
      return cuts[0];  // shortcut if no etaBins
    }
    const auto eta = std::atanh(std::cos(trackState.predicted()[eBoundTheta]));
    const auto abseta = std::abs(eta);
    std::size_t bin = 0;
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

}  // namespace Acts


