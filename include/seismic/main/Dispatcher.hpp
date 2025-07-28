#pragma once

#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"

/**
 * @file Dispatcher.hpp
 * @brief Main simulation dispatcher
 */

namespace seismic {

/**
 * @brief Run complete seismic simulation
 * @param config Simulation configuration
 */
void run_simulation(const SimulationConfig& config);

} // namespace seismic
