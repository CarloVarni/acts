// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <iterator>
#include <memory>

namespace Acts {

class SeedCreator {
 public:
  template <typename spacepoint_iterator_t, typename external_spacepoint_t,
            template <typename...> typename container_t,
            typename platform_t = void*>
  static void createSeeds(
      const Acts::SpacePointGridConfig& gridCfg,
      const Acts::SeedfinderConfig<external_spacepoint_t>& finderCfg,
      spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
      std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
      std::function<std::pair<Acts::Vector3, Acts::Vector2>(
          const external_spacepoint_t&, float, float, float)>
          extractCovariance);
};

template <typename spacepoint_iterator_t, typename external_spacepoint_t,
          template <typename...> typename container_t, typename platform_t>
void SeedCreator::createSeeds(
    const Acts::SpacePointGridConfig& gridCfg,
    const Acts::SeedfinderConfig<external_spacepoint_t>& finderCfg,
    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
    std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
    std::function<std::pair<Acts::Vector3, Acts::Vector2>(
        const external_spacepoint_t&, float, float, float)>
        extractCovariance) {
  auto bottomBinFinder =
      std::make_shared<Acts::BinFinder<external_spacepoint_t>>();
  auto topBinFinder =
      std::make_shared<Acts::BinFinder<external_spacepoint_t>>();
  auto grid =
      Acts::SpacePointGridCreator::createGrid<external_spacepoint_t>(gridCfg);
  auto spacePointsGrouping = Acts::BinnedSPGroup<external_spacepoint_t>(
      spBegin, spEnd, extractCovariance, bottomBinFinder, topBinFinder,
      std::move(grid), finderCfg);

  auto finder = Acts::Seedfinder<external_spacepoint_t, platform_t>(finderCfg);

  static thread_local typename decltype(finder)::State state;

  auto group = spacePointsGrouping.begin();
  auto groupEnd = spacePointsGrouping.end();
  for (; group != groupEnd; ++group) {
    finder.createSeedsForGroup(state, outIt, group.bottom(), group.middle(),
                               group.top());
  }
}

}  // namespace Acts
