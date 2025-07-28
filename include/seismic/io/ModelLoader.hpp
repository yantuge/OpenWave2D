#pragma once

#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"
#include "seismic/core/Grid.hpp"

/**
 * @file ModelLoader.hpp
 * @brief Model loading dispatcher and builders
 */

namespace seismic {

/**
 * @brief Load model parameters into grid based on configuration
 * @param grid Grid object to populate
 * @param config Model configuration
 */
void load_model(Grid& grid, const ModelConfig& config);

/**
 * @brief Build homogeneous model
 * @param grid Grid object to populate
 * @param config Model configuration
 */
void build_homogeneous_model(Grid& grid, const ModelConfig& config);

/**
 * @brief Build anisotropic model
 * @param grid Grid object to populate  
 * @param config Model configuration
 */
void build_anisotropic_model(Grid& grid, const ModelConfig& config);

/**
 * @brief Load model from binary files
 * @param grid Grid object to populate
 * @param config Model configuration
 */
void load_model_from_binary(Grid& grid, const ModelConfig& config);

} // namespace seismic
