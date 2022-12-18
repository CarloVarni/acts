// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CandidatesForSpM.hpp"

namespace Acts {

  CandidatesForSpM::CandidatesForSpM(std::size_t n)
    : m_max_size(n),
      m_SpB(std::numeric_limits<std::size_t>::max())
  {}
  
  void CandidatesForSpM::push(std::size_t SpT, float weight, float zOrigin)
  {
    // if there is still space, add anything
    if (m_storage.size() < m_max_size) {
      addToCollection(m_SpB, SpT, weight, zOrigin);
      return;
    }
    
    // if no space, replace one if quality is enough
    // compare to element with lower weight
    auto& lower_element = m_storage.top();
    if (weight <= std::get<Components::WEIGHT>(lower_element))
      return;
    
    // remove element with lower weight and add this one
    m_storage.pop();
    addToCollection(m_SpB, SpT, weight, zOrigin);
  }
  
  void CandidatesForSpM::addToCollection(std::size_t SpB, std::size_t SpT, float weight, float zOrigin)
  {
    auto toAdd = std::make_tuple(SpB, SpT, weight, zOrigin);
    m_storage.push( toAdd );
  }  
    
} //namespace

