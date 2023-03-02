// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/Holders.hpp" // To be moved somewhere else

#pragma once

namespace Acts::detail {

  // Grid iterator without custom Z-looping
  template <typename grid_type>
  class GridIterator {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = typename grid_type::value_type;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;


    GridIterator() = delete;

    // Constructor
    GridIterator(grid_type& grid, std::size_t index);

    // Never take the ownerships
    GridIterator(GridIterator&&) noexcept = delete;
    
    // Do not allow Copy operations
    GridIterator(const GridIterator&) = delete;
    GridIterator& operator=(const GridIterator&) = delete;

    // Move operations
    GridIterator& operator=(GridIterator&&) noexcept = default;

    // Destructor
    ~GridIterator() = default;

    // Do not allow post-
    GridIterator operator++(int) = delete;
    GridIterator operator--(int) = delete;

    // Only allow pre-
    GridIterator& operator++();
    GridIterator& operator--();

    // Comparisons
    bool operator==(const GridIterator& other) const;
    bool operator!=(const GridIterator& other) const;
    bool operator<(const GridIterator& other) const;
    bool operator>(const GridIterator& other) const;
    bool operator<=(const GridIterator& other) const;
    bool operator>=(const GridIterator& other) const;
    
    typename grid_type::const_reference operator*() const;

  private:
    // The grid
    Acts::detail_tc::RefHolder<grid_type> m_grid;
    // Corresponds to the global index in thr grid
    std::size_t m_index;
  };

  template<typename grid_type>
  GridIterator<grid_type>::GridIterator(grid_type& grid,
					std::size_t index)
    : m_grid( grid ),
      m_index( index )
  {}
  
  template<typename grid_type>
  inline
  GridIterator<grid_type>&
  GridIterator<grid_type>::operator++()
  {
    ++m_index;
    return *this;
  }

  template<typename grid_type>
  inline
  GridIterator<grid_type>&
  GridIterator<grid_type>::operator--()
  {
    --m_index;
    return *this;
  }

  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator==(const GridIterator<grid_type>& other) const
  {
    return m_grid.ptr == other.m_grid.ptr
      and m_index == other.m_index;
  }
  
  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator!=(const GridIterator<grid_type>& other) const
  { return not (*this == other); }

  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator<(const GridIterator<grid_type>& other) const
  { return m_index < other.m_index; }

  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator>(const GridIterator<grid_type>& other) const
  { return not (*this < other)
      and not (*this == other); }

  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator<=(const GridIterator<grid_type>& other) const
  { return not (*this > other); }

  template<typename grid_type>
  inline
  bool
  GridIterator<grid_type>::operator>=(const GridIterator<grid_type>& other) const
  { return not (*this < other); }
  
  template<typename grid_type>
  inline
  typename grid_type::const_reference
  GridIterator<grid_type>::operator*() const
  { return m_grid->at(m_index); }

}
