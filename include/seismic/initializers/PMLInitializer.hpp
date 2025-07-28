#pragma once

#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"
#include "seismic/core/Grid.hpp"

/**
 * @file PMLInitializer.hpp
 * @brief PML profile initialization dispatcher
 */

namespace seismic {

/**
 * @brief Initialize PML profiles based on configuration
 * @param grid Grid object to initialize
 * @param config PML configuration
 * @param model_config Model configuration for velocity information
 * @param dt Time step size for coefficient calculation
 */
void initialize_pml(Grid& grid, const PMLConfig& config, const ModelConfig& model_config, Real dt);

/**
 * @brief Compute CPML damping profiles
 * @param grid Grid object to initialize
 * @param config PML configuration
 * @param quasi_cp_max Maximum P-wave velocity for d0 calculation
 * @param dt Time step size for coefficient calculation
 */
void compute_cpml_profiles(Grid& grid, const PMLConfig& config, Real quasi_cp_max, Real dt);

/**
 * @brief Compute ADE-PML damping profiles
 * @param grid Grid object to initialize
 * @param config PML configuration
 * @param quasi_cp_max Maximum P-wave velocity for d0 calculation
 * @param dt Time step size for coefficient calculation
 */
void compute_adepml_profiles(Grid& grid, const PMLConfig& config, Real quasi_cp_max, Real dt);

} // namespace seismic
