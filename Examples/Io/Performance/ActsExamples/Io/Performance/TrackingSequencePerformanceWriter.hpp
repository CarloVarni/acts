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

#include <mutex>
#include <string>
#include <vector>

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

namespace {
  // Measurement
  using Measurement = ::Acts::BoundVariantMeasurement<ActsExamples::IndexSourceLink>;
  using MeasurementContainer = std::vector<Measurement>;

  // space point
  using SpacePoint = ActsExamples::SimSpacePoint;
  using SpacePointContainer = std::vector<SpacePoint>;

  // seed
  using Seed = Acts::Seed<SpacePoint>;
  using SeedContainer = std::vector<Seed>;

  //  using ParticleContainer = ActsExamples::SimParticleContainer;
  //  using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

  // proto tracks
  using ProtoTrack = ActsExamples::ProtoTrack;
  using ProtoTrackContainer = ActsExamples::ProtoTrackContainer;
}  // namespace 





class TFile;
class TTree;

namespace ActsExamples {
  class TrackingSequencePerformanceWriter final : public WriterT<TrajectoriesContainer> {
  private:
    class seed_info;
    class traj_info;
 public:
  struct Config {
    /// Output Filename
    std::string filePath = "performance_track_seeding.root";
    // output dir
    std::string outputDir;
    /// Tracking geometry for transformation lookup
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Input hit to particles map.
    std::string inputMeasurementParticlesMap;
    /// Input Seed collection
    std::string inputSeeds;
    // Input Space Points
    std::vector<std::string> inputSpacePoints;
    /// Input Proto Tracks
    std::string inputProtoTracks;
    /// Input reconstructed tracks collection.
    std::string inputTracks;
    /// Input truth particles collection.
    std::string inputParticles;

    /// Plot tools
    EffPlotTool::Config effPlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
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


  SeedContainer createSeeds(const ProtoTrackContainer& protoTracks,
			    const SpacePointContainer& spacePoints) const;
  
  SpacePoint localToGlobal(const AlgorithmContext& ctx,
			   const Measurement& meas) const;

  void writetoFile(std::unordered_map<Index, Index>& map_hit_sp,
		   const SpacePointContainer&sps,
		   const std::vector<seed_info>& seedInfoCollection,
		   const std::vector<traj_info>& TrackInfoCollection,
		   std::ofstream& mos);


  private:
  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for efficiency
  /*
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;

  /// Plot tool for duplication rate
  DuplicationPlotTool m_duplicationPlotTool;
  DuplicationPlotTool::DuplicationPlotCache m_duplicationPlotCache{};

  size_t m_nTotalSeeds = 0;
  size_t m_nTotalMatchedSeeds = 0;
  size_t m_nTotalParticles = 0;
  size_t m_nTotalMatchedParticles = 0;
  size_t m_nTotalDuplicatedParticles = 0;
  */

  private:

  class seed_info {
  public:
    seed_info() = delete;
    seed_info(const Seed& seed,
	      Index iseed) 
      : m_measIndexes({
	  seed.sp().at(0)->measurementIndex(), 
	    seed.sp().at(1)->measurementIndex(), 
	    seed.sp().at(2)->measurementIndex()
	    }),
	m_index(iseed) {}
    ~seed_info() = default;
    
    std::string dump() const {
      std::string output = "seed: ";
      output += std::to_string(m_measIndexes.at(0))+ " ";
      output += std::to_string(m_measIndexes.at(1))+ " ";
      output += std::to_string(m_measIndexes.at(2));
      return output;
    }
    
    const Index index() const { return m_index; }
    const std::vector<Index> measIndexes() const { return m_measIndexes; }

  private:
    const std::vector<Index> m_measIndexes;
    const Index m_index;
  };




  class traj_info {
  public:
    traj_info() = delete;
    traj_info(const Trajectories& traj, 
	      unsigned int index)
      : m_index(index) 
    {

      const auto& trackTips = traj.tips();
      const auto& mj = traj.multiTrajectory();

      for (const auto& tip : trackTips) {
	std::vector<Index> indexes;
	std::vector<int> hitType;

	mj.visitBackwards(tip, [&](const auto& state) {
	    if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) return;
	    //and 
	    //		not state.typeFlags().test(Acts::TrackStateFlag::OutlierFlag) )
	    //	      return;
	    
	    if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
	      hitType.push_back(hit_type::MEASUREMENT);
	    } else {
	      hitType.push_back(hit_type::OUTLIER);
	    }
	    
	    Index hitIndex = state.uncalibrated().index();
	    indexes.push_back(hitIndex);
	  });
	
	// Reverse and store
	std::vector<Index> toAdd;
	std::vector<int> toAddType;
	for (unsigned int i(0); i<indexes.size(); i++) {
	  toAdd.push_back( indexes.at(indexes.size() - 1 - i) );
	  toAddType.push_back(hitType.at(hitType.size() - 1 -i));
	}
	m_hits.push_back(toAdd);
      }
    }
    ~traj_info() = default;
    
    std::string dump() const {
      std::string output = "";
      for (long unsigned int i(0); i<m_hits.size(); i++) {
	auto& set = m_hits.at(i);
	for (long unsigned int j(0); j<set.size(); j++) {
	  auto gg = set.at(j);
	  //	  if (m_hitType.at(i).at(j) == hit_type::OUTLIER)
	  //	    output += "*";
	  output += std::to_string(gg);
	  //	  if (m_hitType.at(i).at(j) == hit_type::OUTLIER)
	  //	    output += "*";
	  output += " ";
	}
	output += "\n";
      }
      return output;
    }

    bool isSeedUsed(const seed_info& seed) const {
      const std::vector<Index> seedIndexes = seed.measIndexes();

      for (auto& set : m_hits) {
	std::vector<unsigned int>::size_type index = 0;
	for (auto value : set) {
	  if (value == seedIndexes.at(index))
	    index++;
	  if (index == seedIndexes.size())
	    return true;
	}
      }
      return false;
    }
    

    std::vector<std::vector<Index>> hits() const { return m_hits; }    
    unsigned int index() const { return m_index; }

  private:
    enum hit_type : int {MEASUREMENT = 0, OUTLIER = 1};

    const unsigned int m_index;
    std::vector<std::vector<int>> m_hitType;
    std::vector<std::vector<Index>> m_hits;
  };

};

}  // namespace ActsExamples
