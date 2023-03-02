// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Seeding/Neighborhood.hpp"

#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @c BinnedSPGroupIterator Allows to iterate over all groups of bins
/// a provided BinFinder can generate for each bin of a provided SPGrid
template <typename external_spacepoint_t>
class BinnedSPGroupIterator {
 public:
  BinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>* spgrid,
                        BinFinder<external_spacepoint_t>* botBinFinder,
                        BinFinder<external_spacepoint_t>* tBinFinder,
                        const std::vector<size_t>& bins)
  {
    grid = spgrid;
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    phiZbins = grid->numLocalBins();
    phiIndex = 1;
    zIndex = 1;
    customZorder = &bins;
    // if m_bins vector was not defined, use z bin 1 (zIndex) to start the
    // iterator, otherwise use the first value in m_bins vector
    size_t this_zIndex = bins.empty() ? zIndex : bins.front();
    outputIndex = grid->globalBinFromLocalBins({phiIndex, this_zIndex});
    currentBin = NeighborhoodVector{
        grid->globalBinFromLocalBins({phiIndex, this_zIndex})};
    bottomBinIndices = m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
    topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
  }

  BinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>* spgrid,
                        BinFinder<external_spacepoint_t>* botBinFinder,
                        BinFinder<external_spacepoint_t>* tBinFinder,
                        size_t phiInd, size_t zInd,
                        const std::vector<size_t>& bins = {}) {
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    grid = spgrid;
    phiIndex = phiInd;
    zIndex = zInd;
    phiZbins = grid->numLocalBins();
    customZorder = &bins;
    // if m_bins vector was not defined, use the next z bin (zInd), otherwise
    // use the z bin value stored in m_bins vector for a custom order
    size_t this_zIndex =
        bins.empty()
            ? zIndex
            : (zIndex <= phiZbins[1] ? bins.at(zIndex - 1) : bins.back());
    outputIndex = grid->globalBinFromLocalBins({phiIndex, this_zIndex});
    currentBin =
        NeighborhoodVector(grid->globalBinFromLocalBins({phiInd, this_zIndex}));
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      bottomBinIndices =
          m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
    }
  }

    BinnedSPGroupIterator& operator++() {
    if (zIndex < phiZbins[1]) {
      zIndex++;
    } else {
      zIndex = 1;
      phiIndex++;
    }

    size_t this_zIndex = zIndex;
    if (customZorder and not customZorder->empty()) {
      this_zIndex = customZorder->at(this_zIndex - 1);
    }

    // set current & neighbor bins only if bin indices valid
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      currentBin = NeighborhoodVector{
          grid->globalBinFromLocalBins({phiIndex, this_zIndex})};
      bottomBinIndices =
          m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
      outputIndex++;
      return *this;
    }
    phiIndex = phiZbins[0];
    zIndex = phiZbins[1] + 1;
    return *this;
  }

  bool operator==(const BinnedSPGroupIterator& otherState) {
    return (zIndex == otherState.zIndex && phiIndex == otherState.phiIndex);
  }

  bool operator!=(const BinnedSPGroupIterator& otherState) {
    return !(this->operator==(otherState));
  }

  Neighborhood<external_spacepoint_t> middle() {
    return Neighborhood<external_spacepoint_t>(currentBin, grid);
  }

  Neighborhood<external_spacepoint_t> bottom() {
    return Neighborhood<external_spacepoint_t>(bottomBinIndices, grid);
  }

  Neighborhood<external_spacepoint_t> top() {
    return Neighborhood<external_spacepoint_t>(topBinIndices, grid);
  }

 private:
  // middle spacepoint bin
  NeighborhoodVector currentBin;
  NeighborhoodVector bottomBinIndices;
  NeighborhoodVector topBinIndices;
  const SpacePointGrid<external_spacepoint_t>* grid;
  size_t phiIndex = 1;
  size_t zIndex = 1;
  size_t outputIndex = 0;
  std::array<long unsigned int, 2ul> phiZbins;
  BinFinder<external_spacepoint_t>* m_bottomBinFinder;
  BinFinder<external_spacepoint_t>* m_topBinFinder;
  const std::vector<size_t>* customZorder;
};

/// @c BinnedSPGroup Provides access to begin and end BinnedSPGroupIterator
/// for given BinFinders and SpacePointGrid.
/// Fulfills the range_expression interface.
template <typename external_spacepoint_t>
class BinnedSPGroup {
 public:
  BinnedSPGroup() = delete;

  template <typename spacepoint_iterator_t, typename callable_t>
  BinnedSPGroup<external_spacepoint_t>(
      spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
      callable_t&& toGlobal,
      std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> botBinFinder,
      std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> tBinFinder,
      std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
      Acts::Extent& rRangeSPExtent,
      const SeedFinderConfig<external_spacepoint_t>& _config,
      const SeedFinderOptions& _options);

  size_t size() const { return m_binnedSP->size(); }

  BinnedSPGroupIterator<external_spacepoint_t> begin() {
    return BinnedSPGroupIterator<external_spacepoint_t>(
        m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get(),
        m_bins);
  }

  BinnedSPGroupIterator<external_spacepoint_t> end() {
    auto phiZbins = m_binnedSP->numLocalBins();
    return BinnedSPGroupIterator<external_spacepoint_t>(
        m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get(),
        phiZbins[0], phiZbins[1] + 1, m_bins);
  }

 private:
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<external_spacepoint_t>> m_binnedSP;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<BinFinder<external_spacepoint_t>> m_topBinFinder;
  std::shared_ptr<BinFinder<external_spacepoint_t>> m_bottomBinFinder;

  std::vector<size_t> m_bins;
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"
