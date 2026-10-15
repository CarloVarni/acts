// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/examples/make_magnetic_field.hpp"

// Project include(s)
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"

// Command line option include(s).
#include "traccc/options/clusterization.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/logging.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/throughput.hpp"

// I/O include(s).
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/utils.hpp"

// Performance measurement include(s).
#include "traccc/performance/throughput.hpp"
#include "traccc/performance/timer.hpp"
#include "traccc/performance/timing_info.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// Indicators include(s).
#include "traccc/examples/indicators.hpp"

// System include(s).
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <memory>

namespace traccc {

template <typename FULL_CHAIN_ALG>
int hough_transform_seeding_cuda(std::string_view description, int argc, char* argv[],
                  	 	 vecmem::memory_resource* host_mr) {
  std::unique_ptr<const traccc::Logger> prelogger = traccc::getDefaultLogger(
      "HoughTransformSeedingExample", traccc::Logging::Level::INFO);
 
  // Return gracefully.
  return 0;
}

}  // namespace traccc
