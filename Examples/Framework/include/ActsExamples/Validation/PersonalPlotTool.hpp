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

#include <map>
#include <memory>
#include <string>

namespace ActsExamples {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class PersonalPlotTool {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::map<std::string, PlotHelpers::Binning> varBinning = {
      {"Eta", PlotHelpers::Binning("#eta", 40, -4, 4)},
      {"Phi", PlotHelpers::Binning("#phi", 12, -3.15, 3.15)},
      {"Pt", PlotHelpers::Binning("pT [GeV/c]", 10, 0, 10)},
      {"N_Reco_Times", PlotHelpers::Binning("N Reco Times",10, 0, 10)}
    };
  };

  /// @brief Nested Cache struct
  struct PersonalPlotCache {
    TEfficiency* trackEff_vs_pT_vs_eta;
    TEfficiency* trackEff_vs_pT_vs_phi;
    TEfficiency* trackEff_vs_eta_vs_phi;

    TH2F* muon_reco_pT_vs_eta;
    TH2F* muon_reco_pT_vs_phi;
    TH2F* muon_reco_eta_vs_phi;

    TH1F* muon_reco_pt;
    TH1F* muon_reco_eta;
    TH1F* muon_reco_phi;

    TH2F* muon_truth_pT_vs_eta;
    TH2F* muon_truth_pT_vs_phi;
    TH2F* muon_truth_eta_vs_phi;

    TH1F* muon_truth_pt;
    TH1F* muon_truth_eta;
    TH1F* muon_truth_phi;

    TH1F* n_reco_times;
  };
  
  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param lvl Message level declaration
  PersonalPlotTool(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief book the efficiency plots
  ///
  /// @param effPlotCache the cache for efficiency plots
  void book(PersonalPlotCache& effPlotCache) const;

  /// @brief fill efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  /// @param truthParticle the truth Particle
  /// @param status the reconstruction status
  void fill_truth(PersonalPlotCache& effPlotCache,
		  const ActsFatras::Particle& truthParticle, int ntimes,
		  bool status) const;
  
  void fill_reco(PersonalPlotCache& effPlotCache,
		 float, float, float) const;

  /// @brief write the efficiency plots to file
  ///
  /// @param effPlotCache cache object for efficiency plots
  void write(const PersonalPlotCache& effPlotCache) const;

  /// @brief delete the efficiency plots
  ///
  /// @param effPlotCache cache object for efficiency plots
  void clear(PersonalPlotCache& effPlotCache) const;

 private:
  Config m_cfg;                                  ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
