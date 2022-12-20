// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// System include(s)
#include <algorithm>
#include <cmath>
#include <utility>

// VecMem include(s).
#include "vecmem/containers/vector.hpp"

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/CandidatesForSpM.hpp"

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.hpp"
#include "Acts/Plugins/Sycl/Seeding/SeedFinder.hpp"

namespace Acts::Sycl {
template <typename external_spacepoint_t>
SeedFinder<external_spacepoint_t>::SeedFinder(
    Acts::SeedFinderConfig<external_spacepoint_t> config,
    const Acts::SeedFinderOptions& options,
    const Acts::Sycl::DeviceExperimentCuts& cuts,
    Acts::Sycl::QueueWrapper wrappedQueue, vecmem::memory_resource& resource,
    vecmem::memory_resource* device_resource)
    : m_config(config),
      m_options(options),
      m_deviceCuts(cuts),
      m_wrappedQueue(std::move(wrappedQueue)),
      m_resource(&resource),
      m_device_resource(device_resource) {
  auto seedFilterConfig = m_config.seedFilter->getSeedFilterConfig();

  // init m_deviceConfig
  m_deviceConfig = Acts::Sycl::detail::DeviceSeedFinderConfig{
      m_config.deltaRMin,
      m_config.deltaRMax,
      m_config.cotThetaMax,
      m_config.collisionRegionMin,
      m_config.collisionRegionMax,
      m_config.maxScatteringAngle2,
      m_config.sigmaScattering,
      m_options.minHelixDiameter2,
      m_options.pT2perRadius,
      seedFilterConfig.deltaInvHelixDiameter,
      seedFilterConfig.impactWeightFactor,
      seedFilterConfig.deltaRMin,
      seedFilterConfig.compatSeedWeight,
      m_config.impactMax,
      seedFilterConfig.compatSeedLimit,
  };
}

template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Acts::Seed<external_spacepoint_t>>
SeedFinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  // As a first step, we create Arrays of Structures (AoS)
  // that are easily comprehensible by the GPU. This allows us
  // less memory access operations than with simple (float) arrays.

  // Creating VecMem vectors of the space points, linked to the host/shared
  // memory resource They will be filled and passed to CreateSeedsForGroup().
  vecmem::vector<detail::DeviceSpacePoint> deviceBottomSPs(m_resource);
  vecmem::vector<detail::DeviceSpacePoint> deviceMiddleSPs(m_resource);
  vecmem::vector<detail::DeviceSpacePoint> deviceTopSPs(m_resource);

  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> bottomSPvec;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> middleSPvec;
  std::vector<Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  for (auto SP : bottomSPs) {
    bottomSPvec.push_back(SP);
  }
  deviceBottomSPs.reserve(bottomSPvec.size());
  for (auto SP : bottomSPvec) {
    deviceBottomSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                               SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : middleSPs) {
    middleSPvec.push_back(SP);
  }
  deviceMiddleSPs.reserve(middleSPvec.size());
  for (auto SP : middleSPvec) {
    deviceMiddleSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                               SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : topSPs) {
    topSPvec.push_back(SP);
  }
  deviceTopSPs.reserve(topSPvec.size());
  for (auto SP : topSPvec) {
    deviceTopSPs.push_back({SP->x(), SP->y(), SP->z(), SP->radius(),
                            SP->varianceR(), SP->varianceZ()});
  }

  // std::vector<std::vector<detail::SeedData>> seeds;
  std::vector<std::vector<detail::SeedData>> seeds;

  // Call the SYCL seeding algorithm
  createSeedsForGroupSycl(m_wrappedQueue, *m_resource, m_device_resource,
                          m_deviceConfig, m_deviceCuts, deviceBottomSPs,
                          deviceMiddleSPs, deviceTopSPs, seeds);

  // Iterate through seeds returned by the SYCL algorithm and perform the last
  // step of filtering for fixed middle SP.
  std::vector< typename CandidatesForSpM<InternalSpacePoint<external_spacepoint_t>>::output_type > candidates;

  auto sorting_function = [] (const auto& i1, const auto& i2) -> bool
  {
	const auto& [bottom_l1, medium_l1, top_l1, weight_l1, zOrigin_l1,
                     isQuality_l1] = i1;
        const auto& [bottom_l2, medium_l2, top_l2, weight_l2, zOrigin_l2,
                     isQuality_l2] = i2;

       if (weight_l1 != weight_l2)
          return weight_l1 > weight_l2;

       // This is for the case when the weights from different seeds
        // are same. This makes cpu & cuda results same

        // medium is the same for all candidates
        float sum_medium =
            medium_l1->y() * medium_l1->y() + medium_l1->z() * medium_l1->z();

        float seed1_sum = sum_medium;
        float seed2_sum = sum_medium;

        seed1_sum +=
            bottom_l1->y() * bottom_l1->y() + bottom_l1->z() * bottom_l1->z();
        seed1_sum += top_l1->y() * top_l1->y() + top_l1->z() * top_l1->z();

        seed2_sum +=
            bottom_l2->y() * bottom_l2->y() + bottom_l2->z() * bottom_l2->z();
        seed2_sum += top_l2->y() * top_l2->y() + top_l2->z() * top_l2->z();

        return seed1_sum > seed2_sum;
  };

  for (size_t mi = 0; mi < seeds.size(); ++mi) {
    candidates.clear();
    for (size_t j = 0; j < seeds[mi].size(); ++j) {
      auto& bottomSP = *(bottomSPvec[seeds[mi][j].bottom]);
      auto& middleSP = *(middleSPvec[mi]);
      auto& topSP = *(topSPvec[seeds[mi][j].top]);
      float weight = seeds[mi][j].weight;

      candidates.emplace_back( &bottomSP, &middleSP, topSP&, weight, 0);
    }
    std::sort(candidates.begin(), candidates.end(), sorting_function);
    int numQualitySeeds = 0;  // not used but needs to be fixed
    m_config.seedFilter->filterSeeds_1SpFixed(candidates, numQualitySeeds,
                                              std::back_inserter(outputVec));
  }
  return outputVec;
}
}  // namespace Acts::Sycl
