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
  class NeighborhoodIterator;
  
  template<typename external_spacepoint_t>
  using candidate_collection_t =
    std::vector<const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>*>;
    
  template <typename external_spacepoint_t>
  class Neighborhood {
  public:
    friend NeighborhoodIterator<external_spacepoint_t>;

    Neighborhood(candidate_collection_t<external_spacepoint_t>&& candidates)
      : m_candidates( std::move(candidates) )
    {}
    Neighborhood(candidate_collection_t<external_spacepoint_t>& candidates) = delete;

    NeighborhoodIterator<external_spacepoint_t> begin() const {
      return {*this, 0, 0};
    }

    NeighborhoodIterator<external_spacepoint_t> end() const {
      return {*this, m_candidates->size(), 0};
    }

    NeighborhoodIterator<external_spacepoint_t> begin() {
      return {*this, 0, 0};
    }

    NeighborhoodIterator<external_spacepoint_t> end() {
      return {*this, m_candidates->size(), 0};
    }

  private:
    const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>& at(std::size_t collection) const
    { return *(m_candidates->at(collection)); }

    Acts::InternalSpacePoint<external_spacepoint_t>& at(std::size_t collection, std::size_t element)
    { return *(m_candidates.val.at(collection)->at(element).get()); }
    
  private:
    Acts::detail_tc::ValueHolder< candidate_collection_t<external_spacepoint_t> > m_candidates;
  };

  template<typename external_spacepoint_t>
  class NeighborhoodIterator {
  public:
    NeighborhoodIterator(Neighborhood<external_spacepoint_t>&& neighborhood,
			 std::size_t index_collection,
			 std::size_t index_element) = delete;
    NeighborhoodIterator(Neighborhood<external_spacepoint_t>& neighborhood,
			 std::size_t index_collection,
			 std::size_t index_element)
      : m_neighborhood( neighborhood ),
	m_index_collection( index_collection ),
	m_index_element( index_element )
    {}

    bool operator==(const NeighborhoodIterator& other) const
    {
      return m_neighborhood.ptr == other.m_neighborhood.ptr and
	m_index_collection == other.m_index_collection and
	m_index_element == other.m_index_element;
    }

    bool operator!=(const NeighborhoodIterator& other) const
    { return not (*this == other); }
    
    NeighborhoodIterator& operator++() {
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
    const Neighborhood<external_spacepoint_t>& neighborhood() const
    { return *m_neighborhood.ptr; }
    
  private:
    Acts::detail_tc::RefHolder< Neighborhood<external_spacepoint_t> > m_neighborhood;
    std::size_t m_index_collection {0};
    std::size_t m_index_element {0};
  };

}  // namespace Acts

