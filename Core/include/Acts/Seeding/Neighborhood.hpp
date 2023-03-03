// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
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
#include "Acts/EventData/Holders.hpp"

#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

  template<typename T>
  class SimpleNeighborhoodIterator;
  
  template<typename external_spacepoint_t>
  using candidate_collection_t =
    std::vector<const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>;
    
  template <typename external_spacepoint_t>
  class SimpleNeighborhood {
  public:
    friend SimpleNeighborhoodIterator<external_spacepoint_t>;

    SimpleNeighborhood(candidate_collection_t<external_spacepoint_t>&& candidates) = delete;
    SimpleNeighborhood(candidate_collection_t<external_spacepoint_t>& candidates)
      : m_candidates( candidates )
    {}
    
    SimpleNeighborhoodIterator<external_spacepoint_t> begin() {
      return {*this, 0, 0};
    }

    SimpleNeighborhoodIterator<external_spacepoint_t> end() {
      return {*this, m_candidates->size(), m_candidates->back()->size()};
    }

  private:
    const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>& at(std::size_t collection) const
    { return *(m_candidates->at(collection)); }

    Acts::InternalSpacePoint<external_spacepoint_t>& at(std::size_t collection, std::size_t element)
    { return *(m_candidates.ptr->at(collection)->at(element).get()); }
    
  private:
    Acts::detail_tc::RefHolder< candidate_collection_t<external_spacepoint_t> > m_candidates;
  };

  template<typename external_spacepoint_t>
  class SimpleNeighborhoodIterator {
  public:
    SimpleNeighborhoodIterator(SimpleNeighborhood<external_spacepoint_t>&& neighborhood,
			       std::size_t index_collection,
                               std::size_t index_element) = delete;
    SimpleNeighborhoodIterator(SimpleNeighborhood<external_spacepoint_t>& neighborhood,
			       std::size_t index_collection,
			       std::size_t index_element)
      : m_neighborhood( neighborhood ),
	m_index_collection( index_collection ),
	m_index_element( index_element )
    {}

    bool operator==(const SimpleNeighborhoodIterator& other) const
    {
      return m_neighborhood.ptr == other.m_neighborhood.ptr and
	m_index_collection == other.m_index_collection and
	m_index_element == other.m_index_element;
    }

    bool operator!=(const SimpleNeighborhoodIterator& other) const
    { return not (*this == other); }
    
    SimpleNeighborhoodIterator& operator++() {
      // Increment element
      // If we are at the end of the collection, change collection
      if (++m_index_element == neighborhood().at(m_index_collection).size()) {
	++m_index_collection;
	m_index_element = 0;
      }
      return *this;
    }

    Acts::InternalSpacePoint<external_spacepoint_t>* operator*() {
      return &(m_neighborhood->at(m_index_collection, m_index_element));
    }

  private:
    const SimpleNeighborhood<external_spacepoint_t>& neighborhood() const
    { return *m_neighborhood.ptr; }
    
  private:
    Acts::detail_tc::RefHolder< SimpleNeighborhood<external_spacepoint_t> > m_neighborhood;
    std::size_t m_index_collection;
    std::size_t m_index_element;
  };

  
using NeighborhoodVector = boost::container::small_vector<size_t, 9>;

/// Iterates over the elements of all bins given
/// by the indices parameter in the given SpacePointGrid.
/// Fullfills the forward iterator.
template <typename external_spacepoint_t>
class NeighborhoodIterator {
 public:
  using sp_it_t = typename std::vector<std::unique_ptr<
      InternalSpacePoint<external_spacepoint_t>>>::const_iterator;

  NeighborhoodIterator() = delete;

  // No point in copying the NeighborhoodVector
  NeighborhoodIterator(const NeighborhoodVector& indices,
                       const SpacePointGrid<external_spacepoint_t>* spgrid)
    : m_indices(indices),
      m_curInd(0),
      m_grid(spgrid)
  {
    if (m_indices.size() > m_curInd) {
      m_curIt = std::begin(spgrid->at(m_indices[m_curInd]));
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }

  NeighborhoodIterator(const NeighborhoodVector& indices,
                       const SpacePointGrid<external_spacepoint_t>* spgrid,
                       size_t curInd, sp_it_t curIt)
    : m_curIt(curIt),
      m_indices(indices),
      m_curInd(curInd),
      m_grid(spgrid)
  {
    if (m_indices.size() > m_curInd) {
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }
  
  static NeighborhoodIterator<external_spacepoint_t> begin(
      const NeighborhoodVector& indices,
      const SpacePointGrid<external_spacepoint_t>* spgrid) {
    auto nIt = NeighborhoodIterator<external_spacepoint_t>(indices, spgrid);
    // advance until first non-empty bin or last bin
    if (nIt.m_curIt == nIt.m_binEnd) {
      ++nIt;
    }
    return nIt;
  }

  NeighborhoodIterator(
      const NeighborhoodIterator<external_spacepoint_t>& other)
    : m_curIt( other.m_curIt ),
      m_binEnd( other.m_binEnd ),
      m_indices( other.m_indices ),
      m_curInd( other.m_curInd ),
      m_grid(other.m_grid)
  {}

  void operator++() {
    // if iterator of current Bin not yet at end, increase
    if (m_curIt != m_binEnd) {
      ++m_curIt;
      // return only if end of current bin still not reached
      if (m_curIt != m_binEnd) {
        return;
      }
    }
    // increase bin index m_curInd until you find non-empty bin
    // or until m_curInd >= m_indices.size()-1
    while (m_curIt == m_binEnd and m_indices.size() - 1 > m_curInd) {
      ++m_curInd;
      m_curIt = std::begin(m_grid->at(m_indices[m_curInd]));
      m_binEnd = std::end(m_grid->at(m_indices[m_curInd]));
    }
  }

  InternalSpacePoint<external_spacepoint_t>* operator*() {
    return m_curIt->get();
  }

  bool operator!=(const NeighborhoodIterator<external_spacepoint_t>& other) {
    return m_curIt != other.m_curIt || m_curInd != other.m_curInd;
  }
  
 private:
  // iterators within current bin
  sp_it_t m_curIt;
  sp_it_t m_binEnd;
  // number of bins
  const NeighborhoodVector& m_indices;
  // current bin
  size_t m_curInd;
  const Acts::SpacePointGrid<external_spacepoint_t>* m_grid;
};

/// @c Neighborhood Used to access iterators to access a group of bins
/// returned by a BinFinder.
/// Fulfills the range_expression interface
template <typename external_spacepoint_t>
class Neighborhood {
 public:
  Neighborhood() = delete;
  Neighborhood(NeighborhoodVector&& indices,
	       const SpacePointGrid<external_spacepoint_t>* spgrid) = delete;
  Neighborhood(const NeighborhoodVector& indices,
               const SpacePointGrid<external_spacepoint_t>* spgrid)
    : m_indices(indices),
      m_spgrid(spgrid)
  {}

  NeighborhoodIterator<external_spacepoint_t> begin() const {
    return NeighborhoodIterator<external_spacepoint_t>::begin(m_indices,
                                                              m_spgrid);
  }
  NeighborhoodIterator<external_spacepoint_t> end() const {
    return NeighborhoodIterator<external_spacepoint_t>(
        m_indices, m_spgrid, m_indices.size() - 1,
        std::end(m_spgrid->at(m_indices.back())));
  }

 private:
  const NeighborhoodVector& m_indices;
  const SpacePointGrid<external_spacepoint_t>* m_spgrid;
};

}  // namespace Acts

