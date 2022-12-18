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
  enum Components : int { BSP=0, MSP=1, WEIGHT=2, ZORIGIN=3, QUALITY=4 };
public:
  CandidatesForSpM(const comparator_t& cmp, std::size_t n,
		   std::vector< std::tuple< std::size_t, std::size_t, float, float, bool > >& registry);
  ~CandidatesForSpM() = default;

  void setBottomSp(std::size_t idx);
  void push(std::size_t SpT, float weight, float zOrigin, bool isQuality);
  
private:
  void addToCollection(std::vector< std::tuple< std::size_t, std::size_t, float, float, bool > >& registry,
		       std::size_t SpT, float weight, float zOrigin, bool isQuality);

  void addToCollectionWithIndex(std::vector< std::tuple< std::size_t, std::size_t, float, float, bool > >& registry,
				std::size_t SpT, float weight, float zOrigin, bool isQuality,
				std::size_t index);
  
public:
  std::size_t m_max_size;
  std::size_t m_SpB;
  std::vector< std::tuple< std::size_t, std::size_t, float, float, bool > >& m_registry;
  std::priority_queue<std::size_t, std::vector<std::size_t>, comparator_t> m_storage;
};

}  // namespace Acts
#include "Acts/Seeding/CandidatesForSpM.ipp"
