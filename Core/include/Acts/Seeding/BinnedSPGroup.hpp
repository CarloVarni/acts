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
#include "Acts/EventData/Holders.hpp"

#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {
  template<typename external_spacepoint_t>
  class BinnedSPGroup;
  
  template<typename external_spacepoint_t>
  using candidate_collection_t =
    std::vector<const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>;
  
  /// @c BinnedSPGroupIterator Allows to iterate over all groups of bins
  /// a provided BinFinder can generate for each bin of a provided SPGrid

  /// SpacePointGrid is a very specific structure.
  /// We know it is 2D and what it contains
  /// No need to be too general with this class  
  template <typename external_spacepoint_t>
  class BinnedSPGroupIterator {
  private:
    enum INDEX : int {PHI=0, Z=1};

  public:
    // Never take ownerships
    BinnedSPGroupIterator(BinnedSPGroup<external_spacepoint_t>&& group,
			  std::size_t) = delete;
    BinnedSPGroupIterator(BinnedSPGroup<external_spacepoint_t>& group,
			  std::size_t index)
      : m_group(group),
        m_max_localBins( m_group->m_grid->numLocalBins() )
    {
      if (index == m_group->m_grid->size()) {
	m_current_localBins = m_max_localBins;
      }

      // Go to the next not-empty bin
      findNotEmptyBin();
    }

    BinnedSPGroupIterator(const BinnedSPGroupIterator&) = delete;
    BinnedSPGroupIterator& operator=(const BinnedSPGroupIterator&) = delete;

    BinnedSPGroupIterator(BinnedSPGroupIterator&&) noexcept = default;
    BinnedSPGroupIterator& operator=(BinnedSPGroupIterator&&) noexcept = default;

    ~BinnedSPGroupIterator() = default;

    BinnedSPGroupIterator& operator++()
    {
      // Increase the position by one
      // if we were on the edge, go up one phi bin and reset z bin
      if (++m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z]) {
	++m_current_localBins[INDEX::PHI];
	m_current_localBins[INDEX::Z] = 0;
      }
      
      // Get the next not-empty bin in the grid
      findNotEmptyBin();
      return *this;
    }

    inline bool operator==(const BinnedSPGroupIterator& other) const {
      return m_group.ptr == other.m_group.ptr and
	m_current_localBins[INDEX::PHI] == m_max_localBins[INDEX::PHI] and
	m_current_localBins[INDEX::Z] == m_max_localBins[INDEX::Z];
    }
    inline bool operator!=(const BinnedSPGroupIterator& other) const { return not (*this == other); }

    std::tuple< Neighborhood<external_spacepoint_t>,
		Neighborhood<external_spacepoint_t>,
		Neighborhood<external_spacepoint_t> >
    // std::tuple< candidate_collection_t<external_spacepoint_t>,
    // 		candidate_collection_t<external_spacepoint_t>,
    // 		candidate_collection_t<external_spacepoint_t> >
    operator*() {     
      // Retrieve here - this is the heavy lifting
      // Less expensive then doing it in the operator++
      m_bottomIterators.clear();
      m_middleIterators.clear();
      m_topIterators.clear();
      
      // Middles
      std::size_t global_index = m_group->m_grid->globalBinFromLocalBins({m_current_localBins[INDEX::PHI], m_group->m_bins[m_current_localBins[INDEX::Z]]});
      m_middleIterators.push_back( &(m_group->m_grid->at(global_index)) );
      
      // Bottoms
      auto bottomBinIndices = m_group->m_bottomBinFinder->findBins(m_current_localBins[INDEX::PHI],
								   m_group->m_bins[m_current_localBins[INDEX::Z]],
								   m_group->m_grid.get());
      m_bottomIterators.reserve(bottomBinIndices.size());
      for (auto idx : bottomBinIndices) {
	if (m_group->m_grid->at(idx).size() == 0) continue;
	m_bottomIterators.push_back( &(m_group->m_grid->at(idx)) );
      }
      
      // Tops
      auto topBinIndices = m_group->m_topBinFinder->findBins(m_current_localBins[INDEX::PHI],
							     m_group->m_bins[m_current_localBins[INDEX::Z]],
							     m_group->m_grid.get());
      m_topIterators.reserve(topBinIndices.size());
      for (auto idx : topBinIndices) {
	if (m_group->m_grid->at(idx).size() == 0) continue;
	m_topIterators.push_back( &(m_group->m_grid->at(idx)) );
      }

      //      return std::make_tuple(m_bottomIterators, m_middleIterators, m_topIterators);
      return std::make_tuple(Neighborhood<external_spacepoint_t>(m_bottomIterators),
			     Neighborhood<external_spacepoint_t>(m_middleIterators),
			     Neighborhood<external_spacepoint_t>(m_topIterators));
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
	  
	  std::size_t zBinIndex = m_group->m_bins.size() == 0
	    ? zBin
	    : m_group->m_bins[zBin];
	  std::size_t index = m_group->m_grid->globalBinFromLocalBins({phiBin, zBinIndex});

	  // Check if there are entries in this bin
	  if (m_group->m_grid->at(index).size() == 0) {
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
        
  private:
    /// The group, it contains the grid and the bin finders
    Acts::detail_tc::RefHolder<BinnedSPGroup<external_spacepoint_t>> m_group;
    /// Max Local Bins - limits of the grid
    std::array< std::size_t, 2 > m_max_localBins;
    /// Current Local Bins
    std::array< std::size_t, 2 > m_current_localBins {0, 0};
    // Candidates
    candidate_collection_t<external_spacepoint_t> m_bottomIterators {};
    candidate_collection_t<external_spacepoint_t> m_middleIterators {};
    candidate_collection_t<external_spacepoint_t> m_topIterators {};
  };
  
/// @c BinnedSPGroup Provides access to begin and end BinnedSPGroupIterator
/// for given BinFinders and SpacePointGrid.
/// Fulfills the range_expression interface.
template <typename external_spacepoint_t>
class BinnedSPGroup {
 public:
  friend BinnedSPGroupIterator<external_spacepoint_t>;
  
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

  size_t size() const { return m_grid->size(); }

  BinnedSPGroupIterator<external_spacepoint_t> begin() {
    return {*this, 0};
  }
  
  BinnedSPGroupIterator<external_spacepoint_t> end() {
    return {*this, m_grid->size()};
  }

  // To be removed ?
  const Acts::SpacePointGrid<external_spacepoint_t>* grid() const { return m_grid.get(); }
  
private:
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<external_spacepoint_t>> m_grid;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<const BinFinder<external_spacepoint_t>> m_topBinFinder;
  std::shared_ptr<const BinFinder<external_spacepoint_t>> m_bottomBinFinder;

  std::vector<size_t> m_bins;
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"
