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

  bool sortingCriteria(const Acts::BoundTrackParameters&, const Acts::BoundTrackParameters&) const;
    
private:
  bool hasBeenUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const;
  void setAsUsed(const Acts::MultiTrajectory::TrackStateProxy& trackstate) const;
    
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

inline
bool MikadoSelector::sortingCriteria(const Acts::BoundTrackParameters& parA, const Acts::BoundTrackParameters& parB) const
{ return parA.absoluteMomentum() > parB.absoluteMomentum(); }

} // namespace
