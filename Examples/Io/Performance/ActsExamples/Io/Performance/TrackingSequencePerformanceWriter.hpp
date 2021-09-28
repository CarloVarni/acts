// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Trajectories.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/PersonalPlotTool.hpp"

#include <mutex>
#include <string>
#include <vector>

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"



class TFile;
class TTree;

namespace ActsExamples {

  class TrackingSequencePerformanceWriter final : 
    public WriterT<TrajectoriesContainer> {

  public:
    struct Config {
      /// Output Filename
      std::string fileName = "";
      // output dir
      std::string outputDir = "";
      /// Input reconstructed tracks collection.
      std::string inputMultiTrajectories = "";
      /// Input particles collection.
      std::string inputParticles = "";
      /// Input hit-particles map collection.
      std::string inputMeasurementParticlesMap= "";
      /// Selection criteria
      std::function<bool(float,float,float)> selection_pt_eta_phi = [=] (float, float, float) -> bool { return true; };
      /// Target particle pdg
      std::vector<Acts::PdgParticle> pdg = {};

      EffPlotTool::Config effPlotToolConfig;
      PersonalPlotTool::Config personalPlotToolConfig;
    };
    
    /// Construct from configuration and log level.
    /// @param config The configuration
    /// @param level
    TrackingSequencePerformanceWriter(Config config, 
				      Acts::Logging::Level level);
    ~TrackingSequencePerformanceWriter() override;
    
    /// Finalize plots.
    ProcessCode endRun() final override;
    
  private:
    ProcessCode writeT(const AlgorithmContext& ctx,
		       const TrajectoriesContainer& tracks) final override;
    
  private:
    const Config m_cfg;
    /// Mutex used to protect multi-threaded writes.
    std::mutex m_writeMutex;
    TFile* m_outputFile{nullptr};

    /// Plot tool for efficiency
    EffPlotTool m_effPlotTool;
    EffPlotTool::EffPlotCache m_effPlotCache;

    PersonalPlotTool m_personalPlotTool;
    PersonalPlotTool::PersonalPlotCache m_personalPlotCache;
  };

}  // namespace ActsExamples
