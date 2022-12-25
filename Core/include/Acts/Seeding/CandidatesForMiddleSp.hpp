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
class CandidatesForMiddleSp {
 public:
  using sp_type = external_space_point_t;

  struct candidate {
    candidate(sp_type& b, sp_type& m, sp_type& t,
              float w, float z, bool q) 
      : bottom(&b), middle(&m), top(&t),
	weight(w), zOrigin(z), isQuality(q) {};

    sp_type* bottom;
    sp_type* middle;
    sp_type* top;
    float weight;
    float zOrigin;
    bool isQuality;
  };

  using value_type = candidate;

  CandidatesForMiddleSp();
  ~CandidatesForMiddleSp() = default;

  void setMaxElements(std::size_t n_low, std::size_t n_high);
  void setMiddleSp(sp_type& idx);
  void setBottomSp(sp_type& idx);
  std::vector<value_type> storage() const;

  void push(sp_type& SpT, float weight, float zOrigin, bool isQuality);
  void clear();

  static bool greaterSort(const value_type& i1, const value_type& i2);

 private:
  void push(std::vector<value_type>& storage, std::size_t& n, const std::size_t& n_max,
	    sp_type& SpB, sp_type& SpM, sp_type& SpT, 
	    float weight, float zOrigin, bool isQuality);

  bool exists(const std::size_t& n, const std::size_t& max_size) const;

  void pop(std::vector<value_type>& storage, std::size_t& current_size);
  float top(const std::vector<value_type>& storage) const;
  float weight(const std::vector<value_type>& storage, std::size_t n) const;

  void bubbleup(std::vector<value_type>& storage, std::size_t n);
  void bubbledw(std::vector<value_type>& storage, std::size_t n,
                std::size_t actual_size);

  void addToCollection(std::vector<value_type>& storage, std::size_t& n, 
		       const std::size_t& n_max,
		       value_type&& element);

 private:
  // sizes
  std::size_t m_max_size_high{0};
  std::size_t m_max_size_low{0};
  std::size_t m_n_high{0};
  std::size_t m_n_low{0};

  // space points
  sp_type *m_SpB = nullptr;
  sp_type *m_SpM = nullptr;

  // storage
  // These vectors are sorted as a min heap tree
  // Each node is lower then its childs
  // Thus, it is guaranteed that the lower elements is at the front
  // Sorting criteria is the seed quality

  // storage for candidates with high quality
  std::vector<value_type> m_storage_high {};
  // storage for candidates with low quality
  std::vector<value_type> m_storage_low {};
};

}  // namespace Acts

#include "Acts/Seeding/CandidatesForMiddleSp.ipp"
