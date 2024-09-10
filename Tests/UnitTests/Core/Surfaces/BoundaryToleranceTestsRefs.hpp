// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

// clang-format off
const std::vector<Acts::Vector2> rectVertices = {
    {-2.000000, -1.000000},
    {2.000000, -1.000000},
    {2.000000, 1.000000},
    {-2.000000, 1.000000}
};

const struct {
  double xmin = -2;
  double xmax = 2;
  double ymin = -1;
  double ymax = 1;
} rectDimensions;

const std::vector<Acts::Vector2> rectTestPoints = {
    {-3.00, -2.00}, {-3.00, -1.60}, {-3.00, -1.20}, {-3.00, -0.80}, {-3.00, -0.40},
    {-3.00, 0.00}, {-3.00, 0.40}, {-3.00, 0.80}, {-3.00, 1.20}, {-3.00, 1.60},
    {-3.00, 2.00}, {-2.40, -2.00}, {-2.40, -1.60}, {-2.40, -1.20}, {-2.40, -0.80},
    {-2.40, -0.40}, {-2.40, 0.00}, {-2.40, 0.40}, {-2.40, 0.80}, {-2.40, 1.20},
    {-2.40, 1.60}, {-2.40, 2.00}, {-1.80, -2.00}, {-1.80, -1.60}, {-1.80, -1.20},
    {-1.80, -0.80}, {-1.80, -0.40}, {-1.80, 0.00}, {-1.80, 0.40}, {-1.80, 0.80},
    {-1.80, 1.20}, {-1.80, 1.60}, {-1.80, 2.00}, {-1.20, -2.00}, {-1.20, -1.60},
    {-1.20, -1.20}, {-1.20, -0.80}, {-1.20, -0.40}, {-1.20, 0.00}, {-1.20, 0.40},
    {-1.20, 0.80}, {-1.20, 1.20}, {-1.20, 1.60}, {-1.20, 2.00}, {-0.60, -2.00},
    {-0.60, -1.60}, {-0.60, -1.20}, {-0.60, -0.80}, {-0.60, -0.40}, {-0.60, 0.00},
    {-0.60, 0.40}, {-0.60, 0.80}, {-0.60, 1.20}, {-0.60, 1.60}, {-0.60, 2.00},
    {0.00, -2.00}, {0.00, -1.60}, {0.00, -1.20}, {0.00, -0.80}, {0.00, -0.40},
    {0.00, 0.00}, {0.00, 0.40}, {0.00, 0.80}, {0.00, 1.20}, {0.00, 1.60},
    {0.00, 2.00}, {0.60, -2.00}, {0.60, -1.60}, {0.60, -1.20}, {0.60, -0.80},
    {0.60, -0.40}, {0.60, 0.00}, {0.60, 0.40}, {0.60, 0.80}, {0.60, 1.20},
    {0.60, 1.60}, {0.60, 2.00}, {1.20, -2.00}, {1.20, -1.60}, {1.20, -1.20},
    {1.20, -0.80}, {1.20, -0.40}, {1.20, 0.00}, {1.20, 0.40}, {1.20, 0.80},
    {1.20, 1.20}, {1.20, 1.60}, {1.20, 2.00}, {1.80, -2.00}, {1.80, -1.60},
    {1.80, -1.20}, {1.80, -0.80}, {1.80, -0.40}, {1.80, 0.00}, {1.80, 0.40},
    {1.80, 0.80}, {1.80, 1.20}, {1.80, 1.60}, {1.80, 2.00}, {2.40, -2.00},
    {2.40, -1.60}, {2.40, -1.20}, {2.40, -0.80}, {2.40, -0.40}, {2.40, 0.00},
    {2.40, 0.40}, {2.40, 0.80}, {2.40, 1.20}, {2.40, 1.60}, {2.40, 2.00},
    {3.00, -2.00}, {3.00, -1.60}, {3.00, -1.20}, {3.00, -0.80}, {3.00, -0.40},
    {3.00, 0.00}, {3.00, 0.40}, {3.00, 0.80}, {3.00, 1.20}, {3.00, 1.60},
    {3.00, 2.00}
};
//const std::vector<Acts::Vector2> rectClosestPoints = {
//    {-2.00, -1.00}, {-2.00, -1.00}, {-2.00, -1.00}, {-2.00, -0.80}, {-2.00, -0.40},
//    {-2.00, 0.00}, {-2.00, 0.40}, {-2.00, 0.80}, {-2.00, 1.00}, {-2.00, 1.00},
//    {-2.00, 1.00}, {-2.00, -1.00}, {-2.00, -1.00}, {-2.00, -1.00}, {-2.00, -0.80},
//    {-2.00, -0.40}, {-2.00, 0.00}, {-2.00, 0.40}, {-2.00, 0.80}, {-2.00, 1.00},
//    {-2.00, 1.00}, {-2.00, 1.00}, {-1.80, -1.00}, {-1.80, -1.00}, {-1.80, -1.00},
//    {-2.00, -0.80}, {-2.00, -0.40}, {-2.00, 0.00}, {-2.00, 0.40}, {-1.80, 1.00},
//    {-1.80, 1.00}, {-1.80, 1.00}, {-1.80, 1.00}, {-1.20, -1.00}, {-1.20, -1.00},
//    {-1.20, -1.00}, {-1.20, -1.00}, {-1.20, -1.00}, {-2.00, 0.00}, {-1.20, 1.00},
//    {-1.20, 1.00}, {-1.20, 1.00}, {-1.20, 1.00}, {-1.20, 1.00}, {-0.60, -1.00},
//    {-0.60, -1.00}, {-0.60, -1.00}, {-0.60, -1.00}, {-0.60, -1.00}, {-0.60, -1.00},
//    {-0.60, 1.00}, {-0.60, 1.00}, {-0.60, 1.00}, {-0.60, 1.00}, {-0.60, 1.00},
//    {0.00, -1.00}, {0.00, -1.00}, {0.00, -1.00}, {0.00, -1.00}, {0.00, -1.00},
//    {0.00, -1.00}, {0.00, 1.00}, {0.00, 1.00}, {0.00, 1.00}, {0.00, 1.00},
//    {0.00, 1.00}, {0.60, -1.00}, {0.60, -1.00}, {0.60, -1.00}, {0.60, -1.00},
//    {0.60, -1.00}, {0.60, -1.00}, {0.60, 1.00}, {0.60, 1.00}, {0.60, 1.00},
//    {0.60, 1.00}, {0.60, 1.00}, {1.20, -1.00}, {1.20, -1.00}, {1.20, -1.00},
//    {1.20, -1.00}, {1.20, -1.00}, {2.00, 0.00}, {1.20, 1.00}, {1.20, 1.00},
//    {1.20, 1.00}, {1.20, 1.00}, {1.20, 1.00}, {1.80, -1.00}, {1.80, -1.00},
//    {1.80, -1.00}, {1.80, -1.00}, {2.00, -0.40}, {2.00, 0.00}, {2.00, 0.40},
//    {1.80, 1.00}, {1.80, 1.00}, {1.80, 1.00}, {1.80, 1.00}, {2.00, -1.00},
//    {2.00, -1.00}, {2.00, -1.00}, {2.00, -0.80}, {2.00, -0.40}, {2.00, 0.00},
//    {2.00, 0.40}, {2.00, 0.80}, {2.00, 1.00}, {2.00, 1.00}, {2.00, 1.00},
//    {2.00, -1.00}, {2.00, -1.00}, {2.00, -1.00}, {2.00, -0.80}, {2.00, -0.40},
//    {2.00, 0.00}, {2.00, 0.40}, {2.00, 0.80}, {2.00, 1.00}, {2.00, 1.00},
//    {2.00, 1.00}
//};
const std::vector<double> rectDistances = {
    1.4142135623730951, 1.1661903789690602, 1.019803902718557, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.019803902718557, 1.1661903789690602, 1.4142135623730951,
    1.0770329614269007, 0.7211102550927979, 0.4472135954999578, 0.3999999999999999,
    0.3999999999999999, 0.3999999999999999, 0.3999999999999999, 0.3999999999999999,
    0.4472135954999579, 0.7211102550927979, 1.0770329614269007, 1.0,
    0.6000000000000001, 0.19999999999999996, -0.19999999999999996,
    -0.19999999999999996, -0.19999999999999996, -0.19999999999999996,
    -0.19999999999999973, 0.20000000000000018, 0.6000000000000001, 1.0, 1.0,
    0.6000000000000001, 0.19999999999999996, -0.20000000000000018,
    -0.6000000000000001, -0.7999999999999998, -0.5999999999999996,
    -0.19999999999999973, 0.20000000000000018, 0.6000000000000001, 1.0, 1.0,
    0.6000000000000001, 0.19999999999999996, -0.20000000000000018,
    -0.6000000000000001, -1.0, -0.5999999999999996, -0.19999999999999973,
    0.20000000000000018, 0.6000000000000001, 1.0, 1.0, 0.6000000000000001,
    0.19999999999999996, -0.20000000000000018, -0.6000000000000001, -1.0,
    -0.5999999999999996, -0.19999999999999973, 0.20000000000000018,
    0.6000000000000001, 1.0, 1.0, 0.6000000000000001, 0.19999999999999996,
    -0.20000000000000018, -0.6000000000000001, -1.0, -0.5999999999999996,
    -0.19999999999999973, 0.20000000000000018, 0.6000000000000001, 1.0, 1.0,
    0.6000000000000001, 0.19999999999999996, -0.20000000000000018,
    -0.6000000000000001, -0.7999999999999998, -0.5999999999999996,
    -0.19999999999999973, 0.20000000000000018, 0.6000000000000001, 1.0, 1.0,
    0.6000000000000001, 0.19999999999999996, -0.20000000000000018,
    -0.20000000000000018, -0.20000000000000018, -0.20000000000000018,
    -0.19999999999999973, 0.20000000000000018, 0.6000000000000001, 1.0,
    1.0770329614269007, 0.7211102550927977, 0.44721359549995743,
    0.39999999999999947, 0.39999999999999947, 0.39999999999999947,
    0.39999999999999947, 0.39999999999999947, 0.44721359549995754,
    0.7211102550927977, 1.0770329614269007, 1.4142135623730951, 1.1661903789690602,
    1.019803902718557, 1.0, 1.0, 1.0, 1.0, 1.0, 1.019803902718557,
    1.1661903789690602, 1.4142135623730951
};

const std::vector<Acts::Vector2> rectShiftedVertices = {
    {1.000000, 2.000000},
    {3.000000, 2.000000},
    {3.000000, 4.000000},
    {1.000000, 4.000000}
};

const struct {
  double xmin = 1;
  double xmax = 3;
  double ymin = 2;
  double ymax = 4;
} rectShiftedDimensions;

const std::vector<Acts::Vector2> rectShiftedTestPoints = {
    {0.00, 1.50}, {0.00, 1.80}, {0.00, 2.10}, {0.00, 2.40}, {0.00, 2.70},
    {0.00, 3.00}, {0.00, 3.30}, {0.00, 3.60}, {0.00, 3.90}, {0.00, 4.20},
    {0.00, 4.50}, {0.40, 1.50}, {0.40, 1.80}, {0.40, 2.10}, {0.40, 2.40},
    {0.40, 2.70}, {0.40, 3.00}, {0.40, 3.30}, {0.40, 3.60}, {0.40, 3.90},
    {0.40, 4.20}, {0.40, 4.50}, {0.80, 1.50}, {0.80, 1.80}, {0.80, 2.10},
    {0.80, 2.40}, {0.80, 2.70}, {0.80, 3.00}, {0.80, 3.30}, {0.80, 3.60},
    {0.80, 3.90}, {0.80, 4.20}, {0.80, 4.50}, {1.20, 1.50}, {1.20, 1.80},
    {1.20, 2.10}, {1.20, 2.40}, {1.20, 2.70}, {1.20, 3.00}, {1.20, 3.30},
    {1.20, 3.60}, {1.20, 3.90}, {1.20, 4.20}, {1.20, 4.50}, {1.60, 1.50},
    {1.60, 1.80}, {1.60, 2.10}, {1.60, 2.40}, {1.60, 2.70}, {1.60, 3.00},
    {1.60, 3.30}, {1.60, 3.60}, {1.60, 3.90}, {1.60, 4.20}, {1.60, 4.50},
    {2.00, 1.50}, {2.00, 1.80}, {2.00, 2.10}, {2.00, 2.40}, {2.00, 2.70},
    {2.00, 3.00}, {2.00, 3.30}, {2.00, 3.60}, {2.00, 3.90}, {2.00, 4.20},
    {2.00, 4.50}, {2.40, 1.50}, {2.40, 1.80}, {2.40, 2.10}, {2.40, 2.40},
    {2.40, 2.70}, {2.40, 3.00}, {2.40, 3.30}, {2.40, 3.60}, {2.40, 3.90},
    {2.40, 4.20}, {2.40, 4.50}, {2.80, 1.50}, {2.80, 1.80}, {2.80, 2.10},
    {2.80, 2.40}, {2.80, 2.70}, {2.80, 3.00}, {2.80, 3.30}, {2.80, 3.60},
    {2.80, 3.90}, {2.80, 4.20}, {2.80, 4.50}, {3.20, 1.50}, {3.20, 1.80},
    {3.20, 2.10}, {3.20, 2.40}, {3.20, 2.70}, {3.20, 3.00}, {3.20, 3.30},
    {3.20, 3.60}, {3.20, 3.90}, {3.20, 4.20}, {3.20, 4.50}, {3.60, 1.50},
    {3.60, 1.80}, {3.60, 2.10}, {3.60, 2.40}, {3.60, 2.70}, {3.60, 3.00},
    {3.60, 3.30}, {3.60, 3.60}, {3.60, 3.90}, {3.60, 4.20}, {3.60, 4.50},
    {4.00, 1.50}, {4.00, 1.80}, {4.00, 2.10}, {4.00, 2.40}, {4.00, 2.70},
    {4.00, 3.00}, {4.00, 3.30}, {4.00, 3.60}, {4.00, 3.90}, {4.00, 4.20},
    {4.00, 4.50}
};
//const std::vector<Acts::Vector2> rectShiftedClosestPoints = {
//    {1.00, 2.00}, {1.00, 2.00}, {1.00, 2.10}, {1.00, 2.40}, {1.00, 2.70},
//    {1.00, 3.00}, {1.00, 3.30}, {1.00, 3.60}, {1.00, 3.90}, {1.00, 4.00},
//    {1.00, 4.00}, {1.00, 2.00}, {1.00, 2.00}, {1.00, 2.10}, {1.00, 2.40},
//    {1.00, 2.70}, {1.00, 3.00}, {1.00, 3.30}, {1.00, 3.60}, {1.00, 3.90},
//    {1.00, 4.00}, {1.00, 4.00}, {1.00, 2.00}, {1.00, 2.00}, {1.00, 2.10},
//    {1.00, 2.40}, {1.00, 2.70}, {1.00, 3.00}, {1.00, 3.30}, {1.00, 3.60},
//    {1.00, 3.90}, {1.00, 4.00}, {1.00, 4.00}, {1.20, 2.00}, {1.20, 2.00},
//    {1.20, 2.00}, {1.00, 2.40}, {1.00, 2.70}, {1.00, 3.00}, {1.00, 3.30},
//    {1.00, 3.60}, {1.20, 4.00}, {1.20, 4.00}, {1.20, 4.00}, {1.60, 2.00},
//    {1.60, 2.00}, {1.60, 2.00}, {1.60, 2.00}, {1.00, 2.70}, {1.00, 3.00},
//    {1.00, 3.30}, {1.60, 4.00}, {1.60, 4.00}, {1.60, 4.00}, {1.60, 4.00},
//    {2.00, 2.00}, {2.00, 2.00}, {2.00, 2.00}, {2.00, 2.00}, {2.00, 2.00},
//    {2.00, 2.00}, {2.00, 4.00}, {2.00, 4.00}, {2.00, 4.00}, {2.00, 4.00},
//    {2.00, 4.00}, {2.40, 2.00}, {2.40, 2.00}, {2.40, 2.00}, {2.40, 2.00},
//    {3.00, 2.70}, {3.00, 3.00}, {3.00, 3.30}, {2.40, 4.00}, {2.40, 4.00},
//    {2.40, 4.00}, {2.40, 4.00}, {2.80, 2.00}, {2.80, 2.00}, {2.80, 2.00},
//    {3.00, 2.40}, {3.00, 2.70}, {3.00, 3.00}, {3.00, 3.30}, {3.00, 3.60},
//    {2.80, 4.00}, {2.80, 4.00}, {2.80, 4.00}, {3.00, 2.00}, {3.00, 2.00},
//    {3.00, 2.10}, {3.00, 2.40}, {3.00, 2.70}, {3.00, 3.00}, {3.00, 3.30},
//    {3.00, 3.60}, {3.00, 3.90}, {3.00, 4.00}, {3.00, 4.00}, {3.00, 2.00},
//    {3.00, 2.00}, {3.00, 2.10}, {3.00, 2.40}, {3.00, 2.70}, {3.00, 3.00},
//    {3.00, 3.30}, {3.00, 3.60}, {3.00, 3.90}, {3.00, 4.00}, {3.00, 4.00},
//    {3.00, 2.00}, {3.00, 2.00}, {3.00, 2.10}, {3.00, 2.40}, {3.00, 2.70},
//    {3.00, 3.00}, {3.00, 3.30}, {3.00, 3.60}, {3.00, 3.90}, {3.00, 4.00},
//    {3.00, 4.00}
//};
const std::vector<double> rectShiftedDistances = {
    1.118033988749895, 1.019803902718557, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0198039027185568, 1.118033988749895, 0.7810249675906654, 0.6324555320336759,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6324555320336757, 0.7810249675906654,
    0.5385164807134504, 0.28284271247461895, 0.19999999999999996,
    0.19999999999999996, 0.19999999999999996, 0.19999999999999996,
    0.19999999999999996, 0.19999999999999996, 0.19999999999999996,
    0.28284271247461845, 0.5385164807134504, 0.5, 0.19999999999999996,
    -0.10000000000000009, -0.20000000000000018, -0.20000000000000018,
    -0.20000000000000018, -0.20000000000000018, -0.20000000000000018,
    -0.10000000000000009, 0.1999999999999993, 0.5, 0.5, 0.19999999999999996,
    -0.10000000000000009, -0.3999999999999999, -0.6000000000000001,
    -0.6000000000000001, -0.6000000000000001, -0.3999999999999999,
    -0.10000000000000009, 0.1999999999999993, 0.5, 0.5, 0.19999999999999996,
    -0.10000000000000009, -0.3999999999999999, -0.7000000000000002, -1.0,
    -0.7000000000000002, -0.3999999999999999, -0.10000000000000009,
    0.1999999999999993, 0.5, 0.5, 0.19999999999999996, -0.10000000000000009,
    -0.3999999999999999, -0.5999999999999996, -0.5999999999999996,
    -0.5999999999999996, -0.3999999999999999, -0.10000000000000009,
    0.1999999999999993, 0.5, 0.5, 0.19999999999999996, -0.10000000000000009,
    -0.19999999999999973, -0.19999999999999973, -0.19999999999999973,
    -0.19999999999999973, -0.19999999999999973, -0.10000000000000009,
    0.1999999999999993, 0.5, 0.5385164807134505, 0.28284271247461906,
    0.20000000000000018, 0.20000000000000018, 0.20000000000000018,
    0.20000000000000018, 0.20000000000000018, 0.20000000000000018,
    0.20000000000000018, 0.2828427124746186, 0.5385164807134505, 0.7810249675906655,
    0.6324555320336759, 0.6000000000000001, 0.6000000000000001, 0.6000000000000001,
    0.6000000000000001, 0.6000000000000001, 0.6000000000000001, 0.6000000000000001,
    0.6324555320336757, 0.7810249675906655, 1.118033988749895, 1.019803902718557,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0198039027185568, 1.118033988749895
};

// clang-format on
