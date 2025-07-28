#pragma once

#include <string>
#include <stdexcept>
#include <cmath>
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#include "SourcePolicies.hpp"

/**
 * @file SourceTimeFunction.hpp
 * @brief Factory for source time functions
 */

namespace seismic {

/**
 * @brief Calculate source amplitude at given time
 * @param time Current simulation time
 * @param config Source configuration
 * @return Source amplitude value
 */
inline Real calculate_source_amplitude(Real time, const SourceConfig& config) {
    if (config.type == "ricker") {
        return RickerPolicy::get_value(time, config);
    } else if (config.type == "gaussian") {
        return GaussianPolicy::get_value(time, config);
    } else if (config.type == "first_derivative") {
        return FirstDerivativePolicy::get_value(time, config);
    } else {
        throw std::runtime_error("Unknown source type: " + config.type);
    }
}

/**
 * @brief Calculate source force components
 * @param amplitude Source amplitude 
 * @param config Source configuration
 * @param force_x Output X component of force
 * @param force_y Output Y component of force
 */
inline void calculate_source_forces(Real amplitude, const SourceConfig& config,
                                   Real& force_x, Real& force_y) {
    const Real angle_rad = config.angle * M_PI / 180.0;
    force_x = amplitude * std::sin(angle_rad);
    force_y = amplitude * std::cos(angle_rad);
}

} // namespace seismic
