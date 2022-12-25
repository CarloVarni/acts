// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <memory>
#include <tuple>
#include <vector>

namespace Acts {

template <typename external_space_point_t>
struct Triplet {
  Triplet(external_space_point_t& b, external_space_point_t& m, external_space_point_t& t, 
	    float w, float z, bool q)
      : bottom(&b),
        middle(&m),
        top(&t),
        weight(w),
        zOrigin(z),
        isQuality(q){};

  external_space_point_t* bottom;
  external_space_point_t* middle;
  external_space_point_t* top;
  float weight;
  float zOrigin;
  bool isQuality;
};


template <typename external_space_point_t>
class CandidatesForMiddleSp {
 public:
  using value_type = Triplet<external_space_point_t>;

  CandidatesForMiddleSp();
  ~CandidatesForMiddleSp() = default;

  void setMaxElements(std::size_t n_low, std::size_t n_high);
  std::vector<value_type> storage() const;

  void push(external_space_point_t& SpB, external_space_point_t& SpM, external_space_point_t& SpT, 
	    float weight, float zOrigin, bool isQuality);
  void clear();

  static bool greaterSort(const value_type& i1, const value_type& i2);

 private:
  void push(std::vector<value_type>& storage, std::size_t& n,
            const std::size_t& n_max, external_space_point_t& SpB, external_space_point_t& SpM, external_space_point_t& SpT,
            float weight, float zOrigin, bool isQuality);

  bool exists(const std::size_t& n, const std::size_t& max_size) const;

  void pop(std::vector<value_type>& storage, std::size_t& current_size);
  float weight(const std::vector<value_type>& storage, std::size_t n) const;

  void bubbleup(std::vector<value_type>& storage, std::size_t n);
  void bubbledw(std::vector<value_type>& storage, std::size_t n,
                std::size_t actual_size);

  void addToCollection(std::vector<value_type>& storage, std::size_t& n,
                       const std::size_t& n_max, value_type&& element);

 private:
  // sizes
  std::size_t m_max_size_high{0};
  std::size_t m_max_size_low{0};
  std::size_t m_n_high{0};
  std::size_t m_n_low{0};

  // storage
  // These vectors are sorted as a min heap tree
  // Each node is lower then its childs
  // Thus, it is guaranteed that the lower elements is at the front
  // Sorting criteria is the seed quality

  // storage for candidates with high quality
  std::vector<value_type> m_storage_high{};
  // storage for candidates with low quality
  std::vector<value_type> m_storage_low{};
};

}  // namespace Acts

#include "Acts/Seeding/CandidatesForMiddleSp.ipp"
