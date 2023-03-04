// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

  // Neighborhood
  template <typename external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::Neighborhood(candidate_collection_t<external_spacepoint_t>& candidates)
    : m_candidates( candidates )
  {}

  template <typename external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::Neighborhood(Neighborhood<external_spacepoint_t>&& other) noexcept
    : m_candidates( std::exchange(other.m_candidates.ptr, nullptr) )
  {}
  
  template <typename external_spacepoint_t>
  Neighborhood<external_spacepoint_t>&
  Neighborhood<external_spacepoint_t>::operator=(Neighborhood<external_spacepoint_t>&& other) noexcept
  {
    m_candidates.ptr = std::exchange(other.m_candidates.ptr, nullptr);
  }

  template <typename external_spacepoint_t>
  inline
  NeighborhoodIterator<external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::begin() const
  { return {*this, 0, 0}; }

  template <typename external_spacepoint_t>
  inline
  NeighborhoodIterator<external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::end() const
  { return {*this, m_candidates->size(), 0}; }
  
  template <typename external_spacepoint_t>
  inline
  NeighborhoodIterator<external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::begin()
  { return {*this, 0, 0}; }

  template <typename external_spacepoint_t>
  inline
  NeighborhoodIterator<external_spacepoint_t>
  Neighborhood<external_spacepoint_t>::end()
  { return {*this, m_candidates->size(), 0}; }
  
  template <typename external_spacepoint_t>
  inline
  const std::vector<std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>&
  Neighborhood<external_spacepoint_t>::at(std::size_t collection) const
  { return *(m_candidates->at(collection)); }
  
  template <typename external_spacepoint_t>
  inline
  Acts::InternalSpacePoint<external_spacepoint_t>&
  Neighborhood<external_spacepoint_t>::at(std::size_t collection, std::size_t element)
  { return *(m_candidates.ptr->at(collection)->at(element).get()); }
  

  // NeighborhoodIterator
  template<typename external_spacepoint_t>
  NeighborhoodIterator<external_spacepoint_t>::NeighborhoodIterator(Neighborhood<external_spacepoint_t>& neighborhood,
								    std::size_t index_collection,
								    std::size_t index_element)
    : m_neighborhood( neighborhood ),
      m_index_collection( index_collection ),
      m_index_element( index_element )
  {}
  
  template<typename external_spacepoint_t>
  inline
  bool
  NeighborhoodIterator<external_spacepoint_t>::operator==(const NeighborhoodIterator<external_spacepoint_t>& other) const
  {
    return m_neighborhood.ptr == other.m_neighborhood.ptr and
      m_index_collection == other.m_index_collection and
      m_index_element == other.m_index_element;
  }
  
  template<typename external_spacepoint_t>
  inline
  bool
  NeighborhoodIterator<external_spacepoint_t>::operator!=(const NeighborhoodIterator<external_spacepoint_t>& other) const
  { return not (*this == other); }
  
  template<typename external_spacepoint_t>
  inline
  NeighborhoodIterator<external_spacepoint_t>&
  NeighborhoodIterator<external_spacepoint_t>::operator++() {
    // Increment element
    // If we are at the end of the collection, change collection
    if (++m_index_element == neighborhood().at(m_index_collection).size()) {
      ++m_index_collection;
      m_index_element = 0;
    }
    return *this;
  }
  
  template<typename external_spacepoint_t>
  inline
  Acts::InternalSpacePoint<external_spacepoint_t>*
  NeighborhoodIterator<external_spacepoint_t>::operator*()
  { return &(m_neighborhood->at(m_index_collection, m_index_element)); }
  
  template<typename external_spacepoint_t>
  inline
  const Neighborhood<external_spacepoint_t>&
  NeighborhoodIterator<external_spacepoint_t>::neighborhood() const
  { return *m_neighborhood.ptr; }
  
}  // namespace Acts

