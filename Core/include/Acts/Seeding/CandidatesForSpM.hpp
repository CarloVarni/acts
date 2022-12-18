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

template<typename comparator_t>
class CandidatesForSpM {
private:
  enum Components : int { MSP=0, WEIGHT=1, ZORIGIN=2, QUALITY=3 };
public:
  CandidatesForSpM(const comparator_t& cmp, std::size_t n);
  ~CandidatesForSpM() = default;

  void push(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
	    std::size_t SpT, float weight, float zOrigin, bool isQuality);
  
private:
  void addToCollection(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
		       std::size_t SpT, float weight, float zOrigin, bool isQuality);

  void addToCollectionWithIndex(std::vector< std::tuple< std::size_t, float, float, bool > >& registry,
				std::size_t SpT, float weight, float zOrigin, bool isQuality,
				std::size_t index);
  
private:
  std::size_t m_max_size;
  std::priority_queue<std::size_t, std::vector<std::size_t>, comparator_t> m_storage;
};

}  // namespace Acts
#include "Acts/Seeding/CandidatesForSpM.ipp"
