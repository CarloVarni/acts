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

  class CandidatesForSpM {
  public:
    using value_type = std::tuple<std::size_t, std::size_t, float, float>;
    
    enum Components : int {
      BSP=0,
      TSP=1,
      WEIGHT=2,
      ZORIGIN=3
    };

    CandidatesForSpM(std::size_t n);
    ~CandidatesForSpM() = default;
    
    void setBottomSp(std::size_t idx);
    const std::vector<value_type>& storage();

    void push(std::size_t SpT, float weight, float zOrigin);

  private:
    void pop();
    float top() const;
    bool exists(std::size_t) const;
    float weight(std::size_t) const;

    void bubbleup(std::size_t);
    void bubbledw(std::size_t);
    
    void addToCollection(std::size_t SpB, std::size_t SpT, float weight, float zOrigin);
    void insertToCollection(std::size_t SpB, std::size_t SpT, float weight, float zOrigin);
    
  public:
    std::size_t m_max_size;
    std::size_t m_n;
    std::size_t m_SpB;
    std::vector< value_type > m_storage;
  };

  inline
  const std::vector<typename CandidatesForSpM::value_type>&
  CandidatesForSpM::storage()
  { return m_storage; }

  inline void CandidatesForSpM::setBottomSp(std::size_t idx)
  { m_SpB = idx; }

  inline float CandidatesForSpM::top() const
  { return weight(0); }
  
  inline bool CandidatesForSpM::exists(std::size_t n) const
  { return n < m_n; }

  inline float CandidatesForSpM::weight(std::size_t n) const
  { return std::get<Components::WEIGHT>(m_storage[n]); }
  
}  // namespace Acts

