// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iostream>
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

#include "Acts/EventData/Holders.hpp"

namespace Acts {

  // SpacePointGrid is a very specific structure.
  // We know it is 2D and what it contains
  // No need to be too general with this class
  
  template <typename external_spacepoint_t>
  class SimpleBinnedSPGroupIterator {
    enum INDEX : int {PHI=0, Z=1};
  public:
    // Never take ownerships
    SimpleBinnedSPGroupIterator(SpacePointGrid<external_spacepoint_t>&&,
				const BinFinder<external_spacepoint_t>& bottomBinFinder,
                                const BinFinder<external_spacepoint_t>& topBinFinder,
				std::size_t) = delete;
    SimpleBinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>& grid,
				const BinFinder<external_spacepoint_t>& bottomBinFinder,
				const BinFinder<external_spacepoint_t>& topBinFinder,
				const std::vector<std::size_t>& zCustom,
				std::size_t index)
      : m_grid(grid),
	m_bottomBinFinder( bottomBinFinder ),
	m_topBinFinder( topBinFinder ),
        m_max_localBins( m_grid->numLocalBins() ),
	m_customZ( zCustom )
    {
      if (index == m_grid->size()) {
	m_current_localBins = m_max_localBins;
      }

      // make the default z-index looping
      if (zCustom.size() == 0) {
	m_customZ.reserve(m_max_localBins[INDEX::Z]);
	for (std::size_t i(0); i<m_max_localBins[INDEX::Z]; ++i)
	  m_customZ[i] = i;
      }

      // Go to the next not-empty bin
      findNotEmptyBin();
    }

    SimpleBinnedSPGroupIterator& operator--(int) = delete;
    SimpleBinnedSPGroupIterator& operator--() = delete;
    
    SimpleBinnedSPGroupIterator& operator++(int) = delete;
    SimpleBinnedSPGroupIterator& operator++()
    {
      // Increase the position by one
      ++m_current_localBins[INDEX::Z];
      // if we were on the edge, go up one phi bin and reset z bin
      if (m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z]) {
	++m_current_localBins[INDEX::PHI];
	m_current_localBins[INDEX::Z] = 0;
      }
      
      // Get the next not-empty bin in the grid
      findNotEmptyBin();
      return *this;
    }

    inline bool operator==(const SimpleBinnedSPGroupIterator& other) const {
      return m_grid.ptr == other.m_grid.ptr and
	m_current_localBins[INDEX::PHI] == m_max_localBins[INDEX::PHI] and
	m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z];
    }
    inline bool operator!=(const SimpleBinnedSPGroupIterator& other) const { return not (*this == other); }

    std::tuple<std::vector<
		 const std::vector<
		   std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>,
	       std::vector<
		 const std::vector<
		   std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>,
	       std::vector<
		 const std::vector<
		   std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>>
    operator*() {     
      // Retrieve here 
      // Less expensive then doing it in the operator++
      clear();
      
      // Middles
      std::size_t global_index = m_grid->globalBinFromLocalBins({m_current_localBins[INDEX::PHI], m_customZ[m_current_localBins[INDEX::Z]]});
      m_middleIterators.push_back( &m_grid->at(global_index) );
      
      // Bottoms
      auto bottomBinIndices = m_bottomBinFinder->findBins(m_current_localBins[INDEX::PHI],
							  m_customZ[m_current_localBins[INDEX::Z]],
							  m_grid.ptr);
      m_bottomIterators.reserve(bottomBinIndices.size());
      for (auto idx : bottomBinIndices) {
	m_bottomIterators.push_back(&m_grid->at(idx));
      }
      
      // Tops
      auto topBinIndices = m_topBinFinder->findBins(m_current_localBins[INDEX::PHI],
						    m_customZ[m_current_localBins[INDEX::Z]],
						    m_grid.ptr);
      m_topIterators.reserve(topBinIndices.size());
      for (auto idx : topBinIndices) {
	m_topIterators.push_back(&m_grid->at(idx));
      }
      
      return std::make_tuple(m_middleIterators, m_bottomIterators, m_topIterators);
    }
    
  private:
    inline void findNotEmptyBin() {
      // Iterate on the grid till we find a not-empty bin
      // We start from the current bin configuration and move forward
      for (std::size_t phiBin(m_current_localBins[INDEX::PHI]);
	   phiBin < m_max_localBins[INDEX::PHI];
	   ++phiBin) {
	for (std::size_t zBin(m_current_localBins[INDEX::Z]);
	     zBin < m_max_localBins[INDEX::Z];
	     ++zBin) {
	  
	  std::size_t zBinIndex = m_customZ[zBin];
	  std::size_t index = m_grid->globalBinFromLocalBins({phiBin, zBinIndex});

	  // Check if there are entries in this bin
	  if (m_grid->at(index).size() == 0) {
	    continue;
	  }

	  // Set the new current bins
	  m_current_localBins[INDEX::PHI] = phiBin;
	  m_current_localBins[INDEX::Z] = zBin;
	  return;
	}
	// Reset z-index
	m_current_localBins[INDEX::Z] = 0;
      }
      
      // Could find nothing ... setting this to end()
      m_current_localBins = m_max_localBins;
    }

    void clear() {
      m_middleIterators.clear();
      m_bottomIterators.clear();
      m_topIterators.clear();
    }
        
  private:
    // The grid
    Acts::detail_tc::RefHolder<const SpacePointGrid<external_spacepoint_t>> m_grid;
    // Bottom Bin Finder
    Acts::detail_tc::RefHolder<const BinFinder<external_spacepoint_t>> m_bottomBinFinder;
    // Top Bin Finder
    Acts::detail_tc::RefHolder<const BinFinder<external_spacepoint_t>> m_topBinFinder;
    // Max Local Bins
    std::array< std::size_t, 2 > m_max_localBins;
    // Current Local Bins
    std::array< std::size_t, 2 > m_current_localBins {0, 0};
    // Custom z-navigation
    std::vector<std::size_t> m_customZ {};
    // All iterators
    std::vector<const std::vector<
                  std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*
                > m_middleIterators;
    std::vector<const std::vector<
		  std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*
		> m_topIterators;
    std::vector<const std::vector<
		  std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*
		> m_bottomIterators;
  };
  
  
/// @c BinnedSPGroupIterator Allows to iterate over all groups of bins
/// a provided BinFinder can generate for each bin of a provided SPGrid
template <typename external_spacepoint_t>
class BinnedSPGroupIterator {
 public:
  BinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>* spgrid,
                        const BinFinder<external_spacepoint_t>* botBinFinder,
                        const BinFinder<external_spacepoint_t>* tBinFinder,
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
                        const BinFinder<external_spacepoint_t>* botBinFinder,
                        const BinFinder<external_spacepoint_t>* tBinFinder,
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
  const BinFinder<external_spacepoint_t>* m_bottomBinFinder;
  const BinFinder<external_spacepoint_t>* m_topBinFinder;
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
      std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> botBinFinder,
      std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> tBinFinder,
      std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
      Acts::Extent& rRangeSPExtent,
      const SeedFinderConfig<external_spacepoint_t>& _config,
      const SeedFinderOptions& _options);

  size_t size() const { return m_binnedSP->size(); }

  //  BinnedSPGroupIterator
  SimpleBinnedSPGroupIterator<external_spacepoint_t> begin() {
    // return BinnedSPGroupIterator<external_spacepoint_t>(
    //     m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get(),
    //     m_bins);
    return SimpleBinnedSPGroupIterator(*m_binnedSP.get(),
				       *m_bottomBinFinder.get(),
				       *m_topBinFinder.get(),
				       m_bins,
				       0);
  }

  //  BinnedSPGroupIterator
  SimpleBinnedSPGroupIterator<external_spacepoint_t> end() {
    return SimpleBinnedSPGroupIterator(*m_binnedSP.get(),
                                       *m_bottomBinFinder.get(),
                                       *m_topBinFinder.get(),
                                       m_bins,
                                       m_binnedSP->size());
	
    // auto phiZbins = m_binnedSP->numLocalBins();
    // return BinnedSPGroupIterator<external_spacepoint_t>(
    //     m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get(),
    //     phiZbins[0], phiZbins[1] + 1, m_bins);
  }

  const Acts::SpacePointGrid<external_spacepoint_t>* grid() const { return m_binnedSP.get(); }
  
 private:
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<external_spacepoint_t>> m_binnedSP;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<const BinFinder<external_spacepoint_t>> m_topBinFinder;
  std::shared_ptr<const BinFinder<external_spacepoint_t>> m_bottomBinFinder;

  std::vector<size_t> m_bins;
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"
