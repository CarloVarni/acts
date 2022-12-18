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
#include <queue>
#include <vector>
#include <tuple>

namespace Acts {

  class CandidatesForSpM {
    using stored_collection = std::tuple<std::size_t, std::size_t, float, float>;
    
    enum Components : int {
      BSP=0,
      TSP=1,
      WEIGHT=2,
      ZORIGIN=3
    };
    
    struct candidate_greater {
      bool operator()(const stored_collection& lhs,
		      const stored_collection& rhs)
      { return std::get<Components::WEIGHT>(lhs) > std::get<Components::WEIGHT>(rhs); }
    };
    
  public:
    CandidatesForSpM(std::size_t n);
    ~CandidatesForSpM() = default;
    
    void setBottomSp(std::size_t idx);
    void push(std::size_t SpT, float weight, float zOrigin);

    const std::vector<stored_collection> registry()
    {
      std::vector<stored_collection> output;
      while (m_storage.size() != 0) {
	output.push_back( m_storage.top() );
	m_storage.pop();
      }
      return output;
    }
    
  private:
    void addToCollection(std::size_t SpB, std::size_t SpT, float weight, float zOrigin);

  public:
    std::size_t m_max_size;
    std::size_t m_SpB;
    std::priority_queue<stored_collection,
			std::vector<stored_collection>,
			candidate_greater> m_storage;
  };

  inline void CandidatesForSpM::setBottomSp(std::size_t idx) { m_SpB = idx; }
  
}  // namespace Acts

