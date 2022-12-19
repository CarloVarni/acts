// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CandidatesForSpM.hpp"

namespace Acts {

  CandidatesForSpM::CandidatesForSpM()
    : m_max_size(0),
      m_n(0),
      m_SpB(std::numeric_limits<std::size_t>::max()),
      m_SpM(std::numeric_limits<std::size_t>::max())
  {}
  
  void CandidatesForSpM::push(std::size_t SpT, float weight, float zOrigin)
  {
    // if there is still space, add anything
    if (m_n < m_max_size) {
      addToCollection(m_SpB, m_SpM, SpT, weight, zOrigin);
      return;
    }

    // if no space, replace one if quality is enough
    // compare to element with lower weight
    const auto& lower_weight = top();
    if (weight <= lower_weight)
      return;
    
    // remove element with lower weight and add this one  
    pop();
    insertToCollection(m_SpB, m_SpM, SpT, weight, zOrigin);
  }
  
  void CandidatesForSpM::addToCollection(std::size_t SpB, std::size_t SpM, std::size_t SpT, float weight, float zOrigin)
  {
    auto toAdd = std::make_tuple(SpB, SpM, SpT, weight, zOrigin);
    m_storage.push_back( toAdd );
    std::size_t added_index = m_n++;
    bubbleup(added_index);
  }  

  void CandidatesForSpM::insertToCollection(std::size_t SpB, std::size_t SpM, std::size_t SpT, float weight, float zOrigin)
  {
    auto toAdd = std::make_tuple(SpB, SpM, SpT, weight, zOrigin);
    m_storage[m_n] = toAdd;
    std::size_t added_index = m_n++;
    bubbleup(added_index);
  }
  
  void CandidatesForSpM::bubbledw(std::size_t n)
  {
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    float current = weight(n);
    std::size_t left_child = 2 * n + 1;
    std::size_t right_child = 2 * n + 2;
    
    // no left child, we stop
    if (not exists(left_child)) return;
    // no right child, left wins
    if (not exists(right_child)) {
      if (weight(left_child) < current) {
	std::swap(m_storage[n], m_storage[left_child]);
	return bubbledw(left_child);
      }
    }
    // both childs
    // left is smaller
    if (weight(left_child) < weight(right_child)) {
      if (weight(left_child) < current) {
	std::swap(m_storage[n], m_storage[left_child]);
	return bubbledw(left_child);
      }
    }
    // right is smaller
    if (weight(right_child) < current) {
      std::swap(m_storage[n], m_storage[right_child]);
      return bubbledw(right_child);
    }
    
  }
  
  void CandidatesForSpM::bubbleup(std::size_t n)
  {
    if (n == 0) return;
    
    // parent: (n - 1) / 2;
    std::size_t parent_idx = (n - 1) / 2;
    
    const auto& current = m_storage[n];
    const auto& parent = m_storage[parent_idx];
    if (std::get<Components::WEIGHT>(parent) <= std::get<Components::WEIGHT>(current))
      { return; }
    
    std::swap(m_storage[n], m_storage[parent_idx]);
    bubbleup(parent_idx);
  }
  
  void CandidatesForSpM::pop()
  {
    m_storage[0] = m_storage[m_n - 1];
    --m_n;
    bubbledw(0);
  }
  
} //namespace

