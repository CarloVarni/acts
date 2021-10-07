// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class PhysicsPerformancePlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string tool_name = "";
    std::string particle_name = "";
    std::map<std::string, PlotHelpers::Binning> varBinning = {
      {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
      {"Phi", PlotHelpers::Binning("#phi", 12, -3.15, 3.15)},
      {"Pt", PlotHelpers::Binning("pT [GeV/c]", 10, 0, 10)},
      {"N_Reco_Times", PlotHelpers::Binning("N Reco Times",10, 0, 10)},
      {"N_Contributing_Particles", PlotHelpers::Binning("N Contributing Particles",20, 0, 10)}
    };
  };

  /// @brief Nested Cache struct
  struct PhysicsPerformancePlotCache {
    TEfficiency* trackEff_vs_pT_vs_eta;
    TEfficiency* trackEff_vs_pT_vs_phi;
    TEfficiency* trackEff_vs_eta_vs_phi;

    TH1F* particle_pt;
    TH1F* particle_eta;
    TH1F* particle_phi;

    TH2F* particle_pT_vs_eta;
    TH2F* particle_pT_vs_phi;
    TH2F* particle_eta_vs_phi;

    TH1F* n_contributing_particles;
    TH1F* n_reco_times;
  };
  
  /// Constructor
  PhysicsPerformancePlotTool() = delete;
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  PhysicsPerformancePlotTool(const Config& cfg, 
			     Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param effPlotCache the cache for efficiency plots
  void book(PhysicsPerformancePlotCache& effPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill(PhysicsPerformancePlotCache& effPlotCache,
	    const ActsFatras::Particle& truthParticle,
	    int ntimes, bool status) const;

  void fill(PhysicsPerformancePlotCache& effPlotCache,
            const ActsFatras::Particle& truthParticle) const;

  void fill(PhysicsPerformancePlotCache& effPlotCache,
            const TrackParameters& trackParams,
	    std::size_t nContributingParticles) const;

  void fill(PhysicsPerformancePlotCache& effPlotCache,
	    const TrackParameters& trackParams) const;



  void fill_nContributingParticles(PhysicsPerformancePlotCache& effPlotCache,
				   std::size_t) const;

  void fill(PhysicsPerformancePlotCache& effPlotCache,
	    float, float, float) const;
  
  
  /// @brief write the efficiency plots to file
  ///
  /// @param effPlotCache cache object for efficiency plots
  void write(const PhysicsPerformancePlotCache& effPlotCache) const;

  /// @brief delete the efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  void clear(PhysicsPerformancePlotCache& effPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
