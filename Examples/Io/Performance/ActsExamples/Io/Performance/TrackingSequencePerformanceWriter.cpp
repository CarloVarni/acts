// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TrackingSequencePerformanceWriter.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <stdexcept>
#include <unordered_map>
#include <fstream>

#include <TFile.h>

ActsExamples::TrackingSequencePerformanceWriter::TrackingSequencePerformanceWriter(
    ActsExamples::TrackingSequencePerformanceWriter::Config config,
    Acts::Logging::Level level)
    : WriterT(config.inputTracks, "SeedingPerformanceWriter", level),
      m_cfg(std::move(config)) {
      //      m_effPlotTool(m_cfg.effPlotToolConfig, level),
      //      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, level) {
  /*
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  */
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally

  auto path = m_cfg.filePath;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the plot tools
  //  m_effPlotTool.book(m_effPlotCache);
  //  m_duplicationPlotTool.book(m_duplicationPlotCache);
}

ActsExamples::TrackingSequencePerformanceWriter::~TrackingSequencePerformanceWriter() {
  //  m_effPlotTool.clear(m_effPlotCache);
  //  m_duplicationPlotTool.clear(m_duplicationPlotCache);
  

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::TrackingSequencePerformanceWriter::endRun() {
  if (m_outputFile) {
    //    m_outputFile->cd();
    //    m_effPlotTool.write(m_effPlotCache);
    //    m_duplicationPlotTool.write(m_duplicationPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::TrackingSequencePerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoriesContainer& tracks) {


  // Read truth information collections
  // const auto& particles =
  //   ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  // const auto& hitParticlesMap =
  //   ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);





  return ProcessCode::SUCCESS;

  /*
  // output csv file per event
  // open per-event file
  std::string path =
    perEventFilepath(m_cfg.outputDir, "CKFperformance.csv", ctx.eventNumber);
  std::ofstream mos(path, std::ofstream::out | std::ofstream::trunc);
  if (!mos) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }
  */

  // Read truth information collections
  /*
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  */
  
  /*
  // Get the seeds
  SeedContainer seeds;
  SpacePointContainer spacePoints;
  // These are actually the triplets (seeds) in case we are not running real seeding
  ProtoTrackContainer protoTracks;

  // Read Space Points
  std::unordered_map<Index, Index> map_hit_sp;

  for (const auto& isp : m_cfg.inputSpacePoints) {
    const auto& sps = ctx.eventStore.get<SpacePointContainer>(isp);
    std::copy(sps.begin(), sps.end(), std::back_inserter(spacePoints));
  }

  for (Index isp(0); isp<spacePoints.size(); isp++) {    
    Index hit_index = spacePoints.at(isp).measurementIndex();
    if (map_hit_sp.find(hit_index) != map_hit_sp.end()) {
      ACTS_INFO("HIT index already used by another Space Point!");
    }
    map_hit_sp[hit_index] = isp;
  }


  // Get seed and Proto Tracks
  std::vector<seed_info> seedInfoCollection;

  if (not m_cfg.inputSeeds.empty()) {
    seeds = ctx.eventStore.get<SeedContainer>(m_cfg.inputSeeds);
  } else {
    protoTracks = ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
    seeds = createSeeds(protoTracks, spacePoints);
  }

  for (Index iseed(0); iseed < seeds.size(); iseed++) { 
    seed_info info(seeds.at(iseed), iseed);
    seedInfoCollection.push_back(info);
  }
  



  // Track Info
  std::vector<traj_info> TrackInfoCollection;
  TrackInfoCollection.reserve(tracks.size());

  for (unsigned int ind(0); ind<tracks.size(); ind++) {
    traj_info info(tracks.at(ind), ind);
    TrackInfoCollection.push_back(info);
  }





  writetoFile(map_hit_sp, 
	      spacePoints,
	      seedInfoCollection,
	      TrackInfoCollection,
	      mos);
  
  mos.close();
  return ProcessCode::SUCCESS;
  */
}

void 
ActsExamples::TrackingSequencePerformanceWriter::writetoFile(std::unordered_map<Index, Index>& map_hit_sp,
							     const SpacePointContainer& sps,
							     const std::vector<seed_info>& seedInfoCollection,
							     const std::vector<traj_info>& TrackInfoCollection,
							     std::ofstream& mos) {

  // write csv header
  mos << "track_id,seed_id,"
      << "hit_indexes,"
      << "sp_indexes,"
      << "seed_indexes,"
      << "sp_coordinates";

  mos << '\n';


  for (Index itrack(0); itrack < TrackInfoCollection.size(); itrack++) {
    const auto& traj = TrackInfoCollection.at(itrack);

    for (Index iseed(0); iseed < seedInfoCollection.size(); iseed++) {
      const auto& seed = seedInfoCollection.at(iseed);
      if (not traj.isSeedUsed(seed)) continue;

      mos << itrack << ",";
      mos << iseed << ",";
      
      const auto& all_hits = traj.hits();
      for (const std::vector<Index>& inds : all_hits) {
	  for (const Index& ind : inds)
	    mos << ind << " ";
	  mos << ",";
	  
	  for (const Index& ind : inds) {
	    std::cout<<"writing "<< ind << " -> " << map_hit_sp.at(ind) << std::endl;
	    mos << map_hit_sp.at(ind) << " ";
	  }
	  mos << ",";

	  for (const Index& ind : inds) {
	    Index sp_ind = map_hit_sp.at(ind);
	    mos << sps[sp_ind].x() << ";"
		<< sps[sp_ind].y() << ";"
		<< sps[sp_ind].z() << " ";
	  }
	  mos << '\n';
      }
    }
  }

}







SpacePoint 
ActsExamples::TrackingSequencePerformanceWriter::localToGlobal(const AlgorithmContext& ctx,
							       const Measurement& measurement) const {
  const auto& source_link = std::visit([](const auto& arg) -> const IndexSourceLink& 
				       { return arg.sourceLink(); },
				       measurement);
  auto moduleGeoId = source_link.geometryId();
  const Acts::Surface* surface =
    m_cfg.trackingGeometry->findSurface(moduleGeoId);
  
  // Extract local info
  auto [localPos, localCov] = std::visit(
      [](const auto& meas) {
	auto expander = meas.expander();
	Acts::BoundVector par = expander * meas.parameters();
	Acts::BoundSymMatrix cov =
	expander * meas.covariance() * expander.transpose();
	// extract local position
	Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
	// extract local position covariance.
	Acts::SymMatrix2 lcov =
	cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
	return std::make_pair(lpar, lcov);
      },
      measurement);

  // transform local position to global coordinates
  Acts::Vector3 globalFakeMom(1, 1, 1);
  Acts::Vector3 globalPos =
    surface->localToGlobal(ctx.geoContext, localPos, globalFakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
    surface->referenceFrame(ctx.geoContext, globalPos, globalFakeMom);

  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);
  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // compute Jacobian from local coordinates to rho/z
  Acts::ActsMatrix<2, 2> jac =
            jacXyzToRhoZ *
    rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
  // compute rho/z variance
  Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

  SpacePoint output(globalPos, var[0], var[1], source_link.index());
  return output;
}

SeedContainer
ActsExamples::TrackingSequencePerformanceWriter::createSeeds(const ProtoTrackContainer& protoTracks,
							     const SpacePointContainer& spacePoints) const {
  SeedContainer seeds;
  seeds.reserve(protoTracks.size());
  
  std::unordered_map<Index, const SpacePoint*> spMap;
  
  for (const auto& i : spacePoints) {
    spMap.emplace(i.measurementIndex(), &i);
  }
  
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    if (protoTrack.size() < 3) {
      ACTS_WARNING("Proto track " << itrack << " size is less than 3.");
      continue;
    }
    // Space points on the proto track
    std::vector<const SpacePoint*> spacePointsOnTrack;
    spacePointsOnTrack.reserve(protoTrack.size());
    // Loop over the hit index on the proto track to find the space points
    for (const auto& hitIndex : protoTrack) {
      auto it = spMap.find(hitIndex);
      if (it != spMap.end()) {
        spacePointsOnTrack.push_back(it->second);
      }
    }
    // At least three space points are required
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }
    // Sort the space points
    std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
              [](const SpacePoint* lhs, const SpacePoint* rhs) {
                return std::hypot(lhs->r(), lhs->z()) <
		  std::hypot(rhs->r(), rhs->z());
              });

    // Loop over the found space points to find the seed with maxium deltaR
    // betweent the bottom and top space point
    // @todo add the check of deltaZ
    float deltaRMin = 1. * Acts::UnitConstants::mm;
    /// Maximum deltaR between space points in a seed
    float deltaRMax = 100. * Acts::UnitConstants::mm;

    bool seedFound = false;
    std::array<size_t, 3> bestSPIndices;
    double maxDeltaR = std::numeric_limits<double>::min();
    for (size_t ib = 0; ib < spacePointsOnTrack.size() - 2; ++ib) {
      for (size_t im = ib + 1; im < spacePointsOnTrack.size() - 1; ++im) {
        for (size_t it = im + 1; it < spacePointsOnTrack.size(); ++it) {
          double bmDeltaR = std::abs(spacePointsOnTrack[im]->r() -
                                     spacePointsOnTrack[ib]->r());
          double mtDeltaR = std::abs(spacePointsOnTrack[it]->r() -
                                     spacePointsOnTrack[im]->r());
          if (bmDeltaR >= deltaRMin and bmDeltaR <= deltaRMax and
              mtDeltaR >= deltaRMin and mtDeltaR <= deltaRMax) {
            if ((bmDeltaR + mtDeltaR) > maxDeltaR) {
              maxDeltaR = bmDeltaR + mtDeltaR;
              bestSPIndices = {ib, im, it};
              seedFound = true;
            }
          }
        }
      }
    }

    if (seedFound) {
      seeds.emplace_back(*spacePointsOnTrack[bestSPIndices[0]],
                         *spacePointsOnTrack[bestSPIndices[1]],
                         *spacePointsOnTrack[bestSPIndices[2]],
                         spacePointsOnTrack[bestSPIndices[1]]->z());
    }
  }
  return seeds;
}
