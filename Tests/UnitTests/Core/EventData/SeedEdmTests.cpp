// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Seed.hpp"

#include <vector>

namespace Acts::Test {

struct SpacePoint {};

BOOST_AUTO_TEST_CASE(seed_edm_constructors) {
  std::array<Acts::Test::SpacePoint, 5> storage{};
  Acts::Seed<Acts::Test::SpacePoint, 5> seed(storage[0], storage[1], storage[2],
                                             storage[3], storage[4]);
  const std::array<const Acts::Test::SpacePoint*, 5>& sps = seed.sp();
  for (std::size_t i(0ul); i < 5ul; ++i) {
    BOOST_CHECK_NE(sps[i], nullptr);
    BOOST_CHECK_EQUAL(sps[i], &storage[i]);
  }

  // Copy 1
  Acts::Seed<Acts::Test::SpacePoint, 5> seed_copy1(seed);
  const std::array<const Acts::Test::SpacePoint*, 5>& sps_copy1 =
      seed_copy1.sp();
  for (std::size_t i(0ul); i < 5ul; ++i) {
    BOOST_CHECK_NE(sps_copy1[i], nullptr);
    BOOST_CHECK_EQUAL(sps_copy1[i], sps[i]);
  }

  // Copy 2
  Acts::Seed<Acts::Test::SpacePoint, 5> seed_copy2{seed};
  const std::array<const Acts::Test::SpacePoint*, 5>& sps_copy2 =
      seed_copy2.sp();
  for (std::size_t i(0ul); i < 5ul; ++i) {
    BOOST_CHECK_NE(sps_copy2[i], nullptr);
    BOOST_CHECK_EQUAL(sps_copy2[i], sps[i]);
  }

  // Collection
  std::vector<Acts::Seed<Acts::Test::SpacePoint, 5>> seeds{seed};
  // Copy 1
  std::vector<Acts::Seed<Acts::Test::SpacePoint, 5>> seeds_copy1(seeds);
  BOOST_CHECK_EQUAL(seeds_copy1.size(), seeds.size());
  // Copy 2
  std::vector<Acts::Seed<Acts::Test::SpacePoint, 5>> seeds_copy2{seeds};
  BOOST_CHECK_EQUAL(seeds_copy2.size(), seeds.size());
}

BOOST_AUTO_TEST_CASE(seed_edm_default) {
  std::array<Acts::Test::SpacePoint, 3> storage{};
  Acts::Seed<Acts::Test::SpacePoint> seed(storage[0], storage[1], storage[2]);
  const std::array<const Acts::Test::SpacePoint*, 3>& sps = seed.sp();
  for (std::size_t i(0ul); i < 3ul; ++i) {
    BOOST_CHECK_NE(sps[i], nullptr);
    BOOST_CHECK_EQUAL(sps[i], &storage[i]);
  }

  seed.setVertexZ(-1.2f);
  BOOST_CHECK_EQUAL(seed.z(), -1.2f);

  seed.setQuality(345.23f);
  BOOST_CHECK_EQUAL(seed.seedQuality(), 345.23f);
}

BOOST_AUTO_TEST_CASE(seed_edm_3d) {
  std::array<Acts::Test::SpacePoint, 3> storage{};
  Acts::Seed<Acts::Test::SpacePoint, 3> seed(storage[0], storage[1],
                                             storage[2]);
  const std::array<const Acts::Test::SpacePoint*, 3>& sps = seed.sp();
  for (std::size_t i(0ul); i < 3ul; ++i) {
    BOOST_CHECK_NE(sps[i], nullptr);
    BOOST_CHECK_EQUAL(sps[i], &storage[i]);
  }

  seed.setVertexZ(-1.2f);
  BOOST_CHECK_EQUAL(seed.z(), -1.2f);

  seed.setQuality(345.23f);
  BOOST_CHECK_EQUAL(seed.seedQuality(), 345.23f);
}

BOOST_AUTO_TEST_CASE(seed_edm_4d) {
  std::array<Acts::Test::SpacePoint, 4> storage{};
  Acts::Seed<Acts::Test::SpacePoint, 4> seed(storage[0], storage[1], storage[2],
                                             storage[3]);
  const std::array<const Acts::Test::SpacePoint*, 4>& sps = seed.sp();
  for (std::size_t i(0ul); i < 4ul; ++i) {
    BOOST_CHECK_NE(sps[i], nullptr);
    BOOST_CHECK_EQUAL(sps[i], &storage[i]);
  }

  seed.setVertexZ(-1.2f);
  BOOST_CHECK_EQUAL(seed.z(), -1.2f);

  seed.setQuality(345.23f);
  BOOST_CHECK_EQUAL(seed.seedQuality(), 345.23f);
}

}  // namespace Acts::Test
