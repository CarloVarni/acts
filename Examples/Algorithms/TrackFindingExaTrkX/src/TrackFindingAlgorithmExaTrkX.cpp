// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/ExaTrkX/TorchGraphStoreHook.hpp"
#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <numeric>

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

namespace {

class ExamplesEdmHook : public Acts::ExaTrkXHook {
  double m_targetPT = 0.5_GeV;
  std::size_t m_targetSize = 3;

  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> m_truthGraphHook;
  std::unique_ptr<Acts::TorchTruthGraphMetricsHook> m_targetGraphHook;
  std::unique_ptr<Acts::TorchGraphStoreHook> m_graphStoreHook;

  const Acts::Logger& logger() const { return *m_logger; }

  struct HitInfo {
    std::size_t spacePointIndex;
    std::int32_t hitIndex;
  };

 public:
  ExamplesEdmHook(const SimSpacePointContainer& spacepoints,
                  const IndexMultimap<Index>& measHitMap,
                  const SimHitContainer& simhits,
                  const SimParticleContainer& particles,
                  std::size_t targetMinHits, double targetMinPT,
                  const Acts::Logger& logger)
      : m_targetPT(targetMinPT),
        m_targetSize(targetMinHits),
        m_logger(logger.clone("MetricsHook")) {
    // Associate tracks to graph, collect momentum
    std::unordered_map<ActsFatras::Barcode, std::vector<HitInfo>> tracks;

    for (auto i = 0ul; i < spacepoints.size(); ++i) {
      const auto measId = spacepoints[i]
                              .sourceLinks()[0]
                              .template get<IndexSourceLink>()
                              .index();

      auto [a, b] = measHitMap.equal_range(measId);
      for (auto it = a; it != b; ++it) {
        const auto& hit = *simhits.nth(it->second);

        tracks[hit.particleId()].push_back({i, hit.index()});
      }
    }

    // Collect edges for truth graph and target graph
    std::vector<std::int64_t> truthGraph;
    std::vector<std::int64_t> targetGraph;

    for (auto& [pid, track] : tracks) {
      // Sort by hit index, so the edges are connected correctly
      std::sort(track.begin(), track.end(), [](const auto& a, const auto& b) {
        return a.hitIndex < b.hitIndex;
      });

      auto found = particles.find(pid);
      if (found == particles.end()) {
        ACTS_WARNING("Did not find " << pid << ", skip track");
        continue;
      }

      for (auto i = 0ul; i < track.size() - 1; ++i) {
        truthGraph.push_back(track[i].spacePointIndex);
        truthGraph.push_back(track[i + 1].spacePointIndex);

        if (found->transverseMomentum() > m_targetPT &&
            track.size() >= m_targetSize) {
          targetGraph.push_back(track[i].spacePointIndex);
          targetGraph.push_back(track[i + 1].spacePointIndex);
        }
      }
    }

    m_truthGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        truthGraph, logger.clone());
    m_targetGraphHook = std::make_unique<Acts::TorchTruthGraphMetricsHook>(
        targetGraph, logger.clone());
    m_graphStoreHook = std::make_unique<Acts::TorchGraphStoreHook>();
  }

  ~ExamplesEdmHook() {}

  auto storedGraph() const { return m_graphStoreHook->storedGraph(); }

  void operator()(const std::any& nodes, const std::any& edges,
                  const std::any& weights) const override {
    ACTS_INFO("Metrics for total graph:");
    (*m_truthGraphHook)(nodes, edges, weights);
    ACTS_INFO("Metrics for target graph (pT > "
              << m_targetPT / Acts::UnitConstants::GeV
              << " GeV, nHits >= " << m_targetSize << "):");
    (*m_targetGraphHook)(nodes, edges, weights);
    (*m_graphStoreHook)(nodes, edges, weights);
  }
};

// TODO do we have these function in the repo somewhere?
float theta(float r, float z) {
  return std::atan2(r, z);
}
float eta(float r, float z) {
  return -std::log(std::tan(0.5 * theta(r, z)));
}

}  // namespace

ActsExamples::TrackFindingAlgorithmExaTrkX::TrackFindingAlgorithmExaTrkX(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(config)),
      m_pipeline(m_cfg.graphConstructor, m_cfg.edgeClassifiers,
                 m_cfg.trackBuilder, logger().clone()) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }

  // Sanitizer run with dummy input to detect configuration issues
  // TODO This would be quite helpful I think, but currently it does not work
  // in general because the stages do not expose the number of node features.
  // However, this must be addressed anyway when we also want to allow to
  // configure this more flexible with e.g. cluster information as input. So
  // for now, we disable this.
#if 0
  if( m_cfg.sanitize ) {
  Eigen::VectorXf dummyInput = Eigen::VectorXf::Random(3 * 15);
  std::vector<float> dummyInputVec(dummyInput.data(),
                                   dummyInput.data() + dummyInput.size());
  std::vector<int> spacepointIDs;
  std::iota(spacepointIDs.begin(), spacepointIDs.end(), 0);
  
  runPipeline(dummyInputVec, spacepointIDs);
  }
#endif

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);

  m_inputSimHits.maybeInitialize(m_cfg.inputSimHits);
  m_inputParticles.maybeInitialize(m_cfg.inputParticles);
  m_inputMeasurementMap.maybeInitialize(m_cfg.inputMeasurementSimhitsMap);

  m_outputGraph.maybeInitialize(m_cfg.outputGraph);

  // reserve space for timing
  m_timing.classifierTimes.resize(
      m_cfg.edgeClassifiers.size(),
      decltype(m_timing.classifierTimes)::value_type{0.f});

  // Check if we want cluster features but do not have them
  const static std::array clFeatures = {
      NodeFeature::eClusterX, NodeFeature::eClusterY,  NodeFeature::eCellCount,
      NodeFeature::eCellSum,  NodeFeature::eCluster1R, NodeFeature::eCluster2R};

  auto wantClFeatures = std::any_of(
      m_cfg.nodeFeatures.begin(), m_cfg.nodeFeatures.end(), [&](const auto& f) {
        return std::find(clFeatures.begin(), clFeatures.end(), f) !=
               clFeatures.end();
      });

  if (wantClFeatures && !m_inputClusters.isInitialized()) {
    throw std::invalid_argument("Cluster features requested, but not provided");
  }

  if (m_cfg.nodeFeatures.size() != m_cfg.featureScales.size()) {
    throw std::invalid_argument(
        "Number of features mismatches number of scale parameters.");
  }
}

/// Allow access to features with nice names

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);

  auto hook = std::make_unique<Acts::ExaTrkXHook>();
  if (m_inputSimHits.isInitialized() && m_inputMeasurementMap.isInitialized()) {
    hook = std::make_unique<ExamplesEdmHook>(
        spacepoints, m_inputMeasurementMap(ctx), m_inputSimHits(ctx),
        m_inputParticles(ctx), m_cfg.targetMinHits, m_cfg.targetMinPT,
        logger());
  }

  std::optional<ClusterContainer> clusters;
  if (m_inputClusters.isInitialized()) {
    clusters = m_inputClusters(ctx);
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  const std::size_t numSpacepoints = spacepoints.size();
  const std::size_t numFeatures = m_cfg.nodeFeatures.size();
  ACTS_DEBUG("Received " << numSpacepoints << " spacepoints");
  ACTS_DEBUG("Construct " << numFeatures << " node features");

  std::vector<float> features(numSpacepoints * numFeatures);
  std::vector<int> spacepointIDs;

  spacepointIDs.reserve(spacepoints.size());

  for (auto isp = 0ul; isp < numSpacepoints; ++isp) {
    const auto& sp = spacepoints[isp];

    // For now just take the first index since does require one single index
    // per spacepoint
    // TODO does it work for the module map construction to use only the first
    // sp?
    const auto& sl1 = sp.sourceLinks()[0].template get<IndexSourceLink>();

    // TODO this makes it a bit useless, refactor so we do not need to pass this
    // to the pipeline
    spacepointIDs.push_back(isp);

    // This should be fine, because check in constructor
    Cluster* cl1 = clusters ? &clusters->at(sl1.index()) : nullptr;
    Cluster* cl2 = cl1;

    if (sp.sourceLinks().size() == 2) {
      const auto& sl2 = sp.sourceLinks()[1].template get<IndexSourceLink>();
      cl2 = clusters ? &clusters->at(sl2.index()) : nullptr;
    }

    // I would prefer to use a std::span or boost::span here once available
    float* f = features.data() + isp * numFeatures;

    using NF = NodeFeature;

    for (auto ift = 0ul; ift < numFeatures; ++ift) {
      // clang-format off
      switch(m_cfg.nodeFeatures[ift]) {
        break; case NF::eR:           f[ift] = std::hypot(sp.x(), sp.y());
        break; case NF::ePhi:         f[ift] = std::atan2(sp.y(), sp.x());
        break; case NF::eZ:           f[ift] = sp.z();
        break; case NF::eX:           f[ift] = sp.x();
        break; case NF::eY:           f[ift] = sp.y();
        break; case NF::eEta:         f[ift] = eta(std::hypot(sp.x(), sp.y()), sp.z());
        break; case NF::eClusterX:    f[ift] = cl1->sizeLoc0;
        break; case NF::eClusterY:    f[ift] = cl1->sizeLoc1;
        break; case NF::eCellSum:     f[ift] = cl1->sumActivations();
        break; case NF::eCellCount:   f[ift] = cl1->channels.size();
        break; case NF::eCluster1R:   f[ift] = std::hypot(cl1->globalPosition[Acts::ePos0], cl1->globalPosition[Acts::ePos1]);
        break; case NF::eCluster2R:   f[ift] = std::hypot(cl2->globalPosition[Acts::ePos0], cl2->globalPosition[Acts::ePos1]);
        break; case NF::eCluster1Phi: f[ift] = std::atan2(cl1->globalPosition[Acts::ePos1], cl1->globalPosition[Acts::ePos0]);
        break; case NF::eCluster2Phi: f[ift] = std::atan2(cl2->globalPosition[Acts::ePos1], cl2->globalPosition[Acts::ePos0]);
        break; case NF::eCluster1Z:   f[ift] = cl1->globalPosition[Acts::ePos2];
        break; case NF::eCluster2Z:   f[ift] = cl2->globalPosition[Acts::ePos2];
        break; case NF::eCluster1Eta: f[ift] = eta(std::hypot(cl1->globalPosition[Acts::ePos0], cl1->globalPosition[Acts::ePos1]), cl1->globalPosition[Acts::ePos2]);
        break; case NF::eCluster2Eta: f[ift] = eta(std::hypot(cl2->globalPosition[Acts::ePos0], cl2->globalPosition[Acts::ePos1]), cl2->globalPosition[Acts::ePos2]);
      }
      // clang-format on

      f[ift] /= m_cfg.featureScales[ift];
    }
  }

  // Run the pipeline
  const auto trackCandidates = [&]() {
    std::lock_guard<std::mutex> lock(m_mutex);

    Acts::ExaTrkXTiming timing;
    auto res = m_pipeline.run(features, spacepointIDs, *hook, &timing);

    m_timing.graphBuildingTime(timing.graphBuildingTime.count());

    assert(timing.classifierTimes.size() == m_timing.classifierTimes.size());
    for (auto [aggr, a] :
         Acts::zip(m_timing.classifierTimes, timing.classifierTimes)) {
      aggr(a.count());
    }

    m_timing.trackBuildingTime(timing.trackBuildingTime.count());

    return res;
  }();

  ACTS_DEBUG("Done with pipeline, received " << trackCandidates.size()
                                             << " candidates");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());

  int nShortTracks = 0;

  for (auto& x : trackCandidates) {
    if (m_cfg.filterShortTracks && x.size() < 3) {
      nShortTracks++;
      continue;
    }

    ProtoTrack onetrack;
    onetrack.reserve(x.size());

    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Removed " << nShortTracks << " with less then 3 hits");
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  if (auto dhook = dynamic_cast<ExamplesEdmHook*>(&*hook);
      dhook && m_outputGraph.isInitialized()) {
    auto graph = dhook->storedGraph();
    std::transform(
        graph.first.begin(), graph.first.end(), graph.first.begin(),
        [&](const auto& a) -> std::int64_t { return spacepointIDs.at(a); });
    m_outputGraph(ctx, std::move(graph));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode TrackFindingAlgorithmExaTrkX::finalize() {
  namespace ba = boost::accumulators;

  ACTS_INFO("Exa.TrkX timing info");
  {
    const auto& t = m_timing.graphBuildingTime;
    ACTS_INFO("- graph building: " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }
  for (const auto& t : m_timing.classifierTimes) {
    ACTS_INFO("- classifier:     " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }
  {
    const auto& t = m_timing.trackBuildingTime;
    ACTS_INFO("- track building: " << ba::mean(t) << " +- "
                                   << std::sqrt(ba::variance(t)));
  }

  return {};
}
