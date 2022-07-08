// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/MeasurementSelector.hpp"  
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <functional>

namespace ActsExamples {

  class MikadoSelector :
    public Acts::MeasurementSelector {
public:
  using Config = Acts::MeasurementSelector::Config;

  MikadoSelector() = delete;
  MikadoSelector(Config config) = delete;
  MikadoSelector(Config config, std::size_t n) 
    : Acts::MeasurementSelector(config),
      m_used_meas(n, false)
    {}
  virtual ~MikadoSelector() = default;

  virtual
  Acts::Result<std::pair<std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator,
			 std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator>>
  select(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
         bool& isOutlier, Acts::LoggerWrapper logger) const override;

private:
    bool hasBeenUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const;
    void setAsUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const;

    std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator
    selectBestCandidate(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
			std::function<bool(const Acts::MultiTrajectory::TrackStateProxy&)> selection_function) const;
    
private:
  mutable std::vector<bool> m_used_meas;
};

inline
bool
MikadoSelector::hasBeenUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const
{
  const auto& sourceLink =
    static_cast<const IndexSourceLink&>(trackstate.uncalibrated());
  std::size_t idMeas = sourceLink.index();
  return m_used_meas[idMeas];
}

inline
void
MikadoSelector::setAsUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const
{
  const auto& sourceLink =
    static_cast<const IndexSourceLink&>(trackstate.uncalibrated());
  std::size_t idMeas = sourceLink.index();
  if (m_used_meas[idMeas])
    throw std::invalid_argument("Mikado Selector is considering a measurement that has already been used in the past.");
  m_used_meas[idMeas] = true;
}


// inline for now -> templated in the future
inline
std::vector<Acts::MultiTrajectory::TrackStateProxy>::iterator	
MikadoSelector::selectBestCandidate(std::vector<Acts::MultiTrajectory::TrackStateProxy>& candidates,
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
  
  
} // namespace
