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

    /// Constructors
    Neighborhood(candidate_collection_t<external_spacepoint_t>&& candidates) = delete;
    Neighborhood(candidate_collection_t<external_spacepoint_t>& candidates);
    
    /// Forbid copies
    Neighborhood(const Neighborhood&) = delete;//default;
    Neighborhood& operator=(const Neighborhood&) = delete; //default;

    // /// Move operations
    Neighborhood(Neighborhood&& other) noexcept;
    Neighborhood& operator=(Neighborhood&& other) noexcept;

    /// Destructor
    ~Neighborhood() = default;
        
    NeighborhoodIterator<external_spacepoint_t> begin() const;
    NeighborhoodIterator<external_spacepoint_t> end() const;

    NeighborhoodIterator<external_spacepoint_t> begin();
    NeighborhoodIterator<external_spacepoint_t> end();

  private:
    const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>& at(std::size_t collection) const;
    Acts::InternalSpacePoint<external_spacepoint_t>& at(std::size_t collection, std::size_t element);
    
  private:
    Acts::detail_tc::RefHolder< candidate_collection_t<external_spacepoint_t> > m_candidates;
  };

  

  template<typename external_spacepoint_t>
  class NeighborhoodIterator {
  public:
    NeighborhoodIterator(Neighborhood<external_spacepoint_t>&& neighborhood,
			 std::size_t index_collection,
			 std::size_t index_element) = delete;
    NeighborhoodIterator(Neighborhood<external_spacepoint_t>& neighborhood,
			 std::size_t index_collection,
			 std::size_t index_element);

    bool operator==(const NeighborhoodIterator& other) const;
    bool operator!=(const NeighborhoodIterator& other) const;
    
    NeighborhoodIterator& operator++();

    Acts::InternalSpacePoint<external_spacepoint_t>* operator*();

  private:
    const Neighborhood<external_spacepoint_t>& neighborhood() const;
    
  private:
    Acts::detail_tc::RefHolder< Neighborhood<external_spacepoint_t> > m_neighborhood;
    std::size_t m_index_collection;
    std::size_t m_index_element;
  };
  
}  // namespace Acts

#include "Acts/Seeding/Neighborhood.ipp"


