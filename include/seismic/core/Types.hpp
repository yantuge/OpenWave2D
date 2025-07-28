#pragma once

#include <cstddef>

/**
 * @file Types.hpp
 * @brief Fundamental type definitions for the seismic simulation
 */

namespace seismic {

/// Primary floating-point type for all computations
using Real = double;

/// Integer type for grid indices
using Index = int;

/// Size type for dimensions
using Size = std::size_t;

} // namespace seismic
