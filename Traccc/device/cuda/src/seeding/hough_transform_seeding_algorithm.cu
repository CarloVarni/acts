// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Library include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/hough_transform_seeding_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/device/count_doublets.hpp"
#include "traccc/seeding/device/count_grid_capacities.hpp"
#include "traccc/seeding/device/count_triplets.hpp"
#include "traccc/seeding/device/find_doublets.hpp"
#include "traccc/seeding/device/find_triplets.hpp"
#include "traccc/seeding/device/populate_grid.hpp"
#include "traccc/seeding/device/reduce_triplet_counts.hpp"
#include "traccc/seeding/device/select_seeds.hpp"
#include "traccc/seeding/device/update_triplet_weights.hpp"

namespace traccc::cuda {

hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const seedfinder_config& finder_config,
    const spacepoint_grid_config& grid_config,
    const seedfilter_config& filter_config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger, await_function_type await_func)
    : device::hough_transform_seeding_algorithm(finder_config, grid_config,
                                        filter_config, mr, copy,
                                        std::move(logger)),
      cuda::algorithm_base{str, std::move(await_func)} {}


void triplet_seeding_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::cuda
