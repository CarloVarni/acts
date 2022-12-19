// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iostream>

#include <memory>
#include <vector>
#include <tuple>

namespace Acts {

  template<typename external_space_point_t>
  class CandidatesForSpM {
    using sp_type = external_space_point_t*;
    using value_type = std::tuple<sp_type, sp_type, sp_type, float, float>;
    static constexpr sp_type default_value = nullptr;
    
    enum Components : int {
      BSP=0,
      MSP,
      TSP,
      WEIGHT,
      ZORIGIN
    };
    
  public:
    CandidatesForSpM();
    ~CandidatesForSpM() = default;

    void setMaxElements(std::size_t n);
    void setMediumSp(sp_type);
    void setBottomSp(sp_type);
    const std::vector<value_type>& storage() const;

    void push(sp_type SpT, float weight, float zOrigin);
    void clear();
    
  private:
    void pop();
    float top() const;
    bool exists(std::size_t) const;
    float weight(std::size_t) const;

    void bubbleup(std::size_t);
    void bubbledw(std::size_t);
    
    void addToCollection(sp_type SpB, sp_type SpM, sp_type SpT, float weight, float zOrigin);
    void insertToCollection(sp_type SpB, sp_type SpM, sp_type SpT, float weight, float zOrigin);
    
  public:
    std::size_t m_max_size;
    std::size_t m_n;
    // space points
    sp_type m_SpB;
    sp_type m_SpM;
    // This vector is sorted as a min heap tree
    // Each node is lower then its childs
    // Thus, it is guaranteed that the lower elements is at the front
    // Sorting criteria is the seed quality 
    std::vector< value_type > m_storage;
  };

  template<typename external_space_point_t>
  inline
  const std::vector<typename CandidatesForSpM<external_space_point_t>::value_type>&
  CandidatesForSpM<external_space_point_t>::storage() const
  { return m_storage; }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setMaxElements(std::size_t n)
  {
    if (m_storage.capacity() < n ) m_storage.reserve(n);
    m_max_size = n;
  }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setMediumSp(typename CandidatesForSpM<external_space_point_t>::sp_type idx)
  { m_SpM = idx; }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setBottomSp(typename CandidatesForSpM<external_space_point_t>::sp_type idx)
  { m_SpB = idx; }

  template<typename external_space_point_t>
  inline float CandidatesForSpM<external_space_point_t>::top() const
  { return weight(0); }

  template<typename external_space_point_t>
  inline bool CandidatesForSpM<external_space_point_t>::exists(std::size_t n) const
  { return n < m_n; }

  template<typename external_space_point_t>
  inline float CandidatesForSpM<external_space_point_t>::weight(std::size_t n) const
  { return std::get<Components::WEIGHT>(m_storage[n]); }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::clear()
  {
    // do not clear max size, this is set only once
    m_n = 0;
    // clean bottom so that we understand there is a problem
    // if in cicle this is not manullay set
    m_SpB = default_value;
    // do not clear spm
    m_storage.clear();
  }
  
}  // namespace Acts

#include "Acts/Seeding/CandidatesForSpM.ipp"
