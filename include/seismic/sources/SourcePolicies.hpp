#pragma once

#include <cmath>
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"

/**
 * @file SourcePolicies.hpp
 * @brief Source time function policy classes
 */

namespace seismic {

/**
 * @brief Ricker wavelet source policy
 */
struct RickerPolicy {
    static Real get_value(Real time, const SourceConfig& config) {
        const Real a = M_PI * M_PI * config.f0 * config.f0;
        const Real t_shifted = time - config.t0;
        const Real factor_time = a * t_shifted * t_shifted;
        
        return config.factor * (1.0 - 2.0 * factor_time) * std::exp(-factor_time);
    }
};

/**
 * @brief Gaussian source policy
 */
struct GaussianPolicy {
    static Real get_value(Real time, const SourceConfig& config) {
        const Real a = M_PI * M_PI * config.f0 * config.f0;
        const Real t_shifted = time - config.t0;
        const Real factor_time = a * t_shifted * t_shifted;
        
        return config.factor * std::exp(-factor_time);
    }
};

/**
 * @brief First derivative of Gaussian source policy
 */
struct FirstDerivativePolicy {
    static Real get_value(Real time, const SourceConfig& config) {
        const Real a = M_PI * M_PI * config.f0 * config.f0;
        const Real t_shifted = time - config.t0;
        const Real factor_time = a * t_shifted * t_shifted;
        
        return config.factor * (-2.0 * a * t_shifted) * std::exp(-factor_time);
    }
};

} // namespace seismic
