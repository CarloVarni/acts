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
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

#include <stdexcept>
#include <unordered_map>
#include <fstream>

#include <TFile.h>

ActsExamples::TrackingSequencePerformanceWriter::TrackingSequencePerformanceWriter(ActsExamples::TrackingSequencePerformanceWriter::Config config,
										   Acts::Logging::Level level)
  : WriterT(config.inputMultiTrajectories, "SeedingPerformanceWriter", level),
    m_cfg(std::move(config)) ,
    m_effPlotTool(m_cfg.effPlotToolConfig, level),
    m_truthPerformancePlotTool(m_cfg.truthPerformancePlotConfig, level),
    m_recoPerformancePlotTool(m_cfg.recoPerformancePlotConfig, level),
    m_fakePerformancePlotTool(m_cfg.fakePerformancePlotConfig, level),
    m_unmatchedPerformancePlotTool(m_cfg.unmatchedPerformancePlotConfig, level)
{
  
  if (m_cfg.fileName.empty()) {
    throw std::invalid_argument("Missing output file name");
  }
  
  if (m_cfg.outputDir.empty()) {
    throw std::invalid_argument("Missing output directory name");
  }
  
  if (m_cfg.inputMultiTrajectories.empty()) {
    throw std::invalid_argument("Missing input multi-trajectory collection");
  }

  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection name");
  }

  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input measurement particle map name");
  }

  // Create output root file  
  auto path = m_cfg.outputDir + "/" + m_cfg.fileName;
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_truthPerformancePlotTool.book(m_truthPerformancePlotCache);
  m_recoPerformancePlotTool.book(m_recoPerformancePlotCache);
  m_fakePerformancePlotTool.book(m_fakePerformancePlotCache);
  m_unmatchedPerformancePlotTool.book(m_unmatchedPerformancePlotCache);
  
}

ActsExamples::TrackingSequencePerformanceWriter::~TrackingSequencePerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_truthPerformancePlotTool.clear(m_truthPerformancePlotCache);
  m_recoPerformancePlotTool.clear(m_recoPerformancePlotCache);
  m_fakePerformancePlotTool.clear(m_fakePerformancePlotCache);
  m_unmatchedPerformancePlotTool.clear(m_unmatchedPerformancePlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode 
ActsExamples::TrackingSequencePerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_truthPerformancePlotTool.write(m_truthPerformancePlotCache);
    m_recoPerformancePlotTool.write(m_recoPerformancePlotCache);
    m_fakePerformancePlotTool.write(m_fakePerformancePlotCache);
    m_unmatchedPerformancePlotTool.write(m_unmatchedPerformancePlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode 
ActsExamples::TrackingSequencePerformanceWriter::writeT(const AlgorithmContext& ctx, 
							const TrajectoriesContainer& multi_tracks) {  

  // utilities
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  using RecoTrackInfo = std::pair<size_t, Acts::BoundTrackParameters>;

  // Read truth input collections
  const auto& particles =
    ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
    ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);

  
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);


  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;
  // particle list: barcode -> particle
  std::unordered_map<ActsFatras::Barcode,
		     ::ActsFatras::Particle> storageTruthParticles;
  // particle list: barcode -> n reco times
  std::unordered_map<ActsFatras::Barcode,
		     unsigned long int> storageTruthParticlesNRecoTimes;
  // particles labeled as fakes
  std::unordered_set<ActsFatras::Barcode> storageFakeParticles;


  for (const auto& particle : particles) {
    const ActsFatras::Barcode& barcode = particle.particleId();
    const Acts::PdgParticle& pdg = particle.pdg();

    // Select only target particles
    bool isTargetParticle = false;
    for (const auto& t_pdg : m_cfg.pdg) {
      if (pdg == t_pdg) {
	isTargetParticle = true;
	break;
      }
    }

    if (not isTargetParticle)
      continue;

    // check if already stored
    if (storageTruthParticlesNRecoTimes.find(barcode) != storageTruthParticlesNRecoTimes.end())
      throw std::runtime_error("Truth particle already in memory");

    // position and requirements on pt and eta and phi
    const auto pT = particle.transverseMomentum();
    const auto phi = Acts::VectorHelpers::phi(particle.unitDirection());
    const auto eta = Acts::VectorHelpers::eta(particle.unitDirection());

    if (not m_cfg.selection_pt_eta_phi(pT, eta, phi))
      continue;

    storageTruthParticlesNRecoTimes[barcode] = 0;
    storageTruthParticles[barcode] = particle;
  }
  


  
  // run on input multi-trajectories
  for (std::size_t itraj = 0; itraj < multi_tracks.size(); ++itraj) {
    const auto& traj = multi_tracks[itraj];
    
    if (traj.empty()) {
      ACTS_WARNING("Empty trajectories object " << itraj);
      continue;
    }

    const auto& mj = traj.multiTrajectory();
    const auto& trackTips = traj.tips();


    // run on tips
    for (const auto& trackTip : trackTips) {
      // Check if the reco track has fitted track parameters
      if (not traj.hasTrackParameters(trackTip)) {
        ACTS_WARNING(
            "No fitted track parameters for trajectory with entry index = "
            << trackTip);
        continue;
      }

      const TrackParameters& fittedParameters = traj.trackParameters(trackTip);
      // Requirement on the pT of the track
      const auto& momentum = fittedParameters.momentum();
      const auto pT = Acts::VectorHelpers::perp(momentum);
      const auto eta = Acts::VectorHelpers::eta(momentum);
      const auto phi = Acts::VectorHelpers::phi(momentum);

      if (not m_cfg.selection_pt_eta_phi(pT, eta, phi))
	continue;

      identifyContributingParticles(hitParticlesMap, traj, trackTip,
                                    particleHitCounts);
      if (particleHitCounts.empty()) {
        ACTS_WARNING(
            "No truth particle associated with this trajectory with entry "
            "index = "
            << trackTip);
        continue;
      }

      size_t nContributingParticles = particleHitCounts.size();

      const ActsFatras::Barcode& majorityParticleId =
	particleHitCounts.front().particleId;
      size_t nMajorityHits = particleHitCounts.front().hitCount;
      
      // Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

      // Remove fakes
      if (nMajorityHits * 1. / trajState.nMeasurements < m_cfg.fakeThreshold) {
	m_fakePerformancePlotTool.fill(m_fakePerformancePlotCache, fittedParameters); 
	storageFakeParticles.insert(majorityParticleId);
	continue;
      }

      // check if particle is target
      if (storageTruthParticlesNRecoTimes.find(majorityParticleId) == storageTruthParticlesNRecoTimes.end())
	continue;

      m_recoPerformancePlotTool.fill(m_recoPerformancePlotCache, fittedParameters, nContributingParticles);
      storageTruthParticlesNRecoTimes.at(majorityParticleId)++;
    }
  }
  


  // Count how many particles have been reconstructed now
  long unsigned int counter = 0;
  for (const auto& [barcode, ntimes] : storageTruthParticlesNRecoTimes) {

    // get truth particle
    const auto& truth_particle = storageTruthParticles.at(barcode);
    
    // Fill efficiency plots
    m_effPlotTool.fill(m_effPlotCache, truth_particle, ntimes > 0);
    m_truthPerformancePlotTool.fill(m_truthPerformancePlotCache,
				    truth_particle,
				    ntimes, 
				    ntimes > 0);

    if (ntimes > 0) counter++;
  }

  // check un matched
  for (const auto& particle : particles) {
    const ActsFatras::Barcode& barcode = particle.particleId();
    const Acts::PdgParticle& pdg = particle.pdg();

    // Select only target particles
    bool isTargetParticle = false;
    for (const auto& t_pdg : m_cfg.pdg) {
      if (pdg == t_pdg) {
        isTargetParticle = true;
        break;
      }
    }

    if (not isTargetParticle) 
      continue;

    // check if matched
    if (storageTruthParticlesNRecoTimes.find(barcode) != storageTruthParticlesNRecoTimes.end())
      continue;
    // check if labelled as fake
    if (storageFakeParticles.find(barcode) != storageFakeParticles.end())
      continue;

    // position and requirements on pt and eta and phi
    const auto pT = particle.transverseMomentum();
    const auto phi = Acts::VectorHelpers::phi(particle.unitDirection());
    const auto eta = Acts::VectorHelpers::eta(particle.unitDirection());
    
    if (not m_cfg.selection_pt_eta_phi(pT, eta, phi))
      continue;

    m_unmatchedPerformancePlotTool.fill(m_unmatchedPerformancePlotCache, particle);
  }

  ACTS_INFO("Reconstructed Particles: " << counter << " / " << storageTruthParticlesNRecoTimes.size());

  return ProcessCode::SUCCESS;
}

