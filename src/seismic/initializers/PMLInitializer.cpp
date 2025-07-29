#include "seismic/initializers/PMLInitializer.hpp"
#include <cmath>
#include <algorithm>

namespace seismic {

void initialize_pml(Grid& grid, const PMLConfig& config, const ModelConfig& model_config, Real dt) {
    // Determine maximum P-wave velocity for d0 calculation
    Real quasi_cp_max = 0.0;
    
    if (model_config.type == "homogeneous") {
        quasi_cp_max = model_config.vp;
    } else if (model_config.type == "anisotropic") {
        // For anisotropic media, use sqrt(c22/rho) and sqrt(c11/rho)
        quasi_cp_max = std::max(std::sqrt(model_config.c22 / model_config.density), 
                               std::sqrt(model_config.c11 / model_config.density));
    } else {
        // For binary file models, assume reasonable default
        quasi_cp_max = 3000.0; // Default reasonable velocity
    }
    
    if (config.type == "cpml") {
        compute_cpml_profiles(grid, config, quasi_cp_max, dt);
    } else if (config.type == "adepml") {
        compute_adepml_profiles(grid, config, quasi_cp_max, dt);
    }
}

void compute_cpml_profiles(Grid& grid, const PMLConfig& config, Real quasi_cp_max, Real dt) {
    const Real thickness_pml_x = config.npoints_pml * grid.deltax();
    const Real thickness_pml_y = config.npoints_pml * grid.deltay();
    
    // Compute d0 from INRIA report section 6.1
    const Real d0_x = -(config.npower + 1) * quasi_cp_max * std::log(config.rcoef) / (2.0 * thickness_pml_x);
    const Real d0_y = -(config.npower + 1) * quasi_cp_max * std::log(config.rcoef) / (2.0 * thickness_pml_y);
    
    const Real xorigin_left = thickness_pml_x;
    const Real xorigin_right = (grid.nx() - 1) * grid.deltax() - thickness_pml_x;
    const Real yorigin_bottom = thickness_pml_y;
    const Real yorigin_top = (grid.ny() - 1) * grid.deltay() - thickness_pml_y;
    
    // Initialize X-direction profiles
    for (Index i = 0; i < grid.nx(); ++i) {
        const Real xval = grid.deltax() * static_cast<Real>(i);
        
        // Left edge
        if (config.use_pml_xmin) {
            const Real abscissa_in_pml = xorigin_left - xval;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_x;
                
                grid.d_x(i) = d0_x * std::pow(abscissa_normalized, config.npower);
                grid.k_x(i) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                grid.alpha_x(i) = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = xorigin_left - (xval + grid.deltax() / 2.0);
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_x;
                
                grid.d_x_half(i) = d0_x * std::pow(abscissa_normalized_half, config.npower);
                grid.k_x_half(i) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                grid.alpha_x_half(i) = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Right edge
        if (config.use_pml_xmax) {
            const Real abscissa_in_pml = xval - xorigin_right;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_x;
                
                grid.d_x(i) = d0_x * std::pow(abscissa_normalized, config.npower);
                grid.k_x(i) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                grid.alpha_x(i) = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = xval + grid.deltax() / 2.0 - xorigin_right;
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_x;
                
                grid.d_x_half(i) = d0_x * std::pow(abscissa_normalized_half, config.npower);
                grid.k_x_half(i) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                grid.alpha_x_half(i) = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Ensure alpha is non-negative
        grid.alpha_x(i) = std::max(0.0, grid.alpha_x(i));
        grid.alpha_x_half(i) = std::max(0.0, grid.alpha_x_half(i));
        
        // Compute b and a coefficients
        grid.b_x(i) = std::exp(-(grid.d_x(i) / grid.k_x(i) + grid.alpha_x(i)) * dt);
        grid.b_x_half(i) = std::exp(-(grid.d_x_half(i) / grid.k_x_half(i) + grid.alpha_x_half(i)) * dt);
        
        if (std::abs(grid.d_x(i)) > 1.0e-6) {
            grid.a_x(i) = grid.d_x(i) * (grid.b_x(i) - 1.0) / 
                         (grid.k_x(i) * (grid.d_x(i) + grid.k_x(i) * grid.alpha_x(i)));
        }
        if (std::abs(grid.d_x_half(i)) > 1.0e-6) {
            grid.a_x_half(i) = grid.d_x_half(i) * (grid.b_x_half(i) - 1.0) / 
                              (grid.k_x_half(i) * (grid.d_x_half(i) + grid.k_x_half(i) * grid.alpha_x_half(i)));
        }
    }
    
    // Initialize Y-direction profiles
    for (Index j = 0; j < grid.ny(); ++j) {
        const Real yval = grid.deltay() * static_cast<Real>(j);
        
        // Bottom edge
        if (config.use_pml_ymin) {
            const Real abscissa_in_pml = yorigin_bottom - yval;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_y;
                
                grid.d_y(j) = d0_y * std::pow(abscissa_normalized, config.npower);
                grid.k_y(j) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                grid.alpha_y(j) = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = yorigin_bottom - (yval + grid.deltay() / 2.0);
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_y;
                
                grid.d_y_half(j) = d0_y * std::pow(abscissa_normalized_half, config.npower);
                grid.k_y_half(j) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                grid.alpha_y_half(j) = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Top edge
        if (config.use_pml_ymax) {
            const Real abscissa_in_pml = yval - yorigin_top;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_y;
                
                grid.d_y(j) = d0_y * std::pow(abscissa_normalized, config.npower);
                grid.k_y(j) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                grid.alpha_y(j) = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = yval + grid.deltay() / 2.0 - yorigin_top;
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_y;
                
                grid.d_y_half(j) = d0_y * std::pow(abscissa_normalized_half, config.npower);
                grid.k_y_half(j) = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                grid.alpha_y_half(j) = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Ensure alpha is non-negative
        grid.alpha_y(j) = std::max(0.0, grid.alpha_y(j));
        grid.alpha_y_half(j) = std::max(0.0, grid.alpha_y_half(j));
        
        // Compute b and a coefficients
        grid.b_y(j) = std::exp(-(grid.d_y(j) / grid.k_y(j) + grid.alpha_y(j)) * dt);
        grid.b_y_half(j) = std::exp(-(grid.d_y_half(j) / grid.k_y_half(j) + grid.alpha_y_half(j)) * dt);
        
        if (std::abs(grid.d_y(j)) > 1.0e-6) {
            grid.a_y(j) = grid.d_y(j) * (grid.b_y(j) - 1.0) / 
                         (grid.k_y(j) * (grid.d_y(j) + grid.k_y(j) * grid.alpha_y(j)));
        }
        if (std::abs(grid.d_y_half(j)) > 1.0e-6) {
            grid.a_y_half(j) = grid.d_y_half(j) * (grid.b_y_half(j) - 1.0) / 
                              (grid.k_y_half(j) * (grid.d_y_half(j) + grid.k_y_half(j) * grid.alpha_y_half(j)));
        }
    }
}

void compute_adepml_profiles(Grid& grid, const PMLConfig& config, Real quasi_cp_max, Real dt) {
    // First compute basic PML profiles (d, k, alpha) similar to CPML
    const Real thickness_pml_x = config.npoints_pml * grid.deltax();
    const Real thickness_pml_y = config.npoints_pml * grid.deltay();
    
    // Compute d0 from INRIA report section 6.1
    const Real d0_x = -(config.npower + 1) * quasi_cp_max * std::log(config.rcoef) / (2.0 * thickness_pml_x);
    const Real d0_y = -(config.npower + 1) * quasi_cp_max * std::log(config.rcoef) / (2.0 * thickness_pml_y);
    
    const Real xorigin_left = thickness_pml_x;
    const Real xorigin_right = (grid.nx() - 1) * grid.deltax() - thickness_pml_x;
    const Real yorigin_bottom = thickness_pml_y;
    const Real yorigin_top = (grid.ny() - 1) * grid.deltay() - thickness_pml_y;
    
    // RK4 coefficients matching the original Fortran implementation
    const Real rk41[4] = {0.0, 0.5, 0.5, 1.0};
    
    // ADE-PML parameters from configuration (matching original Fortran code)
    const Real epsn = config.epsn;   // coefficient for ADE
    const Real epsn1 = config.epsn1; // coefficient for ADE
    
    // Extended grid dimensions for ADE-PML (including ghost points)
    const Index nx_extended = grid.nx() + 8;  // -4 to nx+4
    const Index ny_extended = grid.ny() + 8;  // -4 to ny+4
    
    // Initialize X-direction profiles for ADE-PML
    for (Index i_ext = 0; i_ext < nx_extended; ++i_ext) {
        // Convert extended index to regular grid coordinate
        const Index i = i_ext - 4;  // Convert from (-4 to nx+4) to (0 to nx+8)
        const Real xval = grid.deltax() * static_cast<Real>(i);
        
        Real d_x_val = 0.0, k_x_val = 1.0, alpha_x_val = 0.0;
        Real d_x_half_val = 0.0, k_x_half_val = 1.0, alpha_x_half_val = 0.0;
        
        // Left edge PML
        if (config.use_pml_xmin && i >= -4 && i <= config.npoints_pml) {
            const Real abscissa_in_pml = xorigin_left - xval;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_x;
                d_x_val = d0_x * std::pow(abscissa_normalized, config.npower);
                k_x_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                alpha_x_val = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = xorigin_left - (xval + grid.deltax() / 2.0);
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_x;
                d_x_half_val = d0_x * std::pow(abscissa_normalized_half, config.npower);
                k_x_half_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                alpha_x_half_val = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Right edge PML  
        if (config.use_pml_xmax && i >= grid.nx() - config.npoints_pml - 1 && i <= grid.nx() + 4) {
            const Real abscissa_in_pml = xval - xorigin_right;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_x;
                d_x_val = d0_x * std::pow(abscissa_normalized, config.npower);
                k_x_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                alpha_x_val = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = xval + grid.deltax() / 2.0 - xorigin_right;
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_x;
                d_x_half_val = d0_x * std::pow(abscissa_normalized_half, config.npower);
                k_x_half_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                alpha_x_half_val = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Ensure alpha is non-negative
        alpha_x_val = std::max(0.0, alpha_x_val);
        alpha_x_half_val = std::max(0.0, alpha_x_half_val);
        
        // Compute ADE-PML coefficients for all 4 RK4 substeps
        for (int rk_substep = 0; rk_substep < 4; ++rk_substep) {
            // Calculate b coefficients (matching original Fortran formulation)
            const Real factor_x = epsn * dt * rk41[rk_substep] * (d_x_val / k_x_val + alpha_x_val);
            const Real factor_x_half = epsn * dt * rk41[rk_substep] * (d_x_half_val / k_x_half_val + alpha_x_half_val);
            const Real factor1_x = epsn1 * dt * rk41[rk_substep] * (d_x_val / k_x_val + alpha_x_val);
            const Real factor1_x_half = epsn1 * dt * rk41[rk_substep] * (d_x_half_val / k_x_half_val + alpha_x_half_val);
            
            grid.b_x_rk4(rk_substep, i_ext) = (1.0 - factor_x) / (1.0 + factor1_x);
            grid.b_x_half_rk4(rk_substep, i_ext) = (1.0 - factor_x_half) / (1.0 + factor1_x_half);
            
            // Calculate a coefficients
            if (std::abs(d_x_val) > 1.0e-6) {
                grid.a_x_rk4(rk_substep, i_ext) = -dt * rk41[rk_substep] * d_x_val / (k_x_val * k_x_val) / (1.0 + factor1_x);
            } else {
                grid.a_x_rk4(rk_substep, i_ext) = 0.0;
            }
            
            if (std::abs(d_x_half_val) > 1.0e-6) {
                grid.a_x_half_rk4(rk_substep, i_ext) = -dt * rk41[rk_substep] * d_x_half_val / (k_x_half_val * k_x_half_val) / (1.0 + factor1_x_half);
            } else {
                grid.a_x_half_rk4(rk_substep, i_ext) = 0.0;
            }
        }
    }
    
    // Initialize Y-direction profiles for ADE-PML
    for (Index j_ext = 0; j_ext < ny_extended; ++j_ext) {
        // Convert extended index to regular grid coordinate
        const Index j = j_ext - 4;  // Convert from (-4 to ny+4) to (0 to ny+8)
        const Real yval = grid.deltay() * static_cast<Real>(j);
        
        Real d_y_val = 0.0, k_y_val = 1.0, alpha_y_val = 0.0;
        Real d_y_half_val = 0.0, k_y_half_val = 1.0, alpha_y_half_val = 0.0;
        
        // Bottom edge PML
        if (config.use_pml_ymin && j >= -4 && j <= config.npoints_pml) {
            const Real abscissa_in_pml = yorigin_bottom - yval;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_y;
                d_y_val = d0_y * std::pow(abscissa_normalized, config.npower);
                k_y_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                alpha_y_val = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = yorigin_bottom - (yval + grid.deltay() / 2.0);
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_y;
                d_y_half_val = d0_y * std::pow(abscissa_normalized_half, config.npower);
                k_y_half_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                alpha_y_half_val = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Top edge PML
        if (config.use_pml_ymax && j >= grid.ny() - config.npoints_pml - 1 && j <= grid.ny() + 4) {
            const Real abscissa_in_pml = yval - yorigin_top;
            if (abscissa_in_pml >= 0.0) {
                const Real abscissa_normalized = abscissa_in_pml / thickness_pml_y;
                d_y_val = d0_y * std::pow(abscissa_normalized, config.npower);
                k_y_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized, config.npower);
                alpha_y_val = config.alpha_max_pml * (1.0 - abscissa_normalized);
            }
            
            const Real abscissa_in_pml_half = yval + grid.deltay() / 2.0 - yorigin_top;
            if (abscissa_in_pml_half >= 0.0) {
                const Real abscissa_normalized_half = abscissa_in_pml_half / thickness_pml_y;
                d_y_half_val = d0_y * std::pow(abscissa_normalized_half, config.npower);
                k_y_half_val = 1.0 + (config.k_max_pml - 1.0) * std::pow(abscissa_normalized_half, config.npower);
                alpha_y_half_val = config.alpha_max_pml * (1.0 - abscissa_normalized_half);
            }
        }
        
        // Ensure alpha is non-negative
        alpha_y_val = std::max(0.0, alpha_y_val);
        alpha_y_half_val = std::max(0.0, alpha_y_half_val);
        
        // Compute ADE-PML coefficients for all 4 RK4 substeps
        for (int rk_substep = 0; rk_substep < 4; ++rk_substep) {
            // Calculate b coefficients
            const Real factor_y = epsn * dt * rk41[rk_substep] * (d_y_val / k_y_val + alpha_y_val);
            const Real factor_y_half = epsn * dt * rk41[rk_substep] * (d_y_half_val / k_y_half_val + alpha_y_half_val);
            const Real factor1_y = epsn1 * dt * rk41[rk_substep] * (d_y_val / k_y_val + alpha_y_val);
            const Real factor1_y_half = epsn1 * dt * rk41[rk_substep] * (d_y_half_val / k_y_half_val + alpha_y_half_val);
            
            grid.b_y_rk4(rk_substep, j_ext) = (1.0 - factor_y) / (1.0 + factor1_y);
            grid.b_y_half_rk4(rk_substep, j_ext) = (1.0 - factor_y_half) / (1.0 + factor1_y_half);
            
            // Calculate a coefficients
            if (std::abs(d_y_val) > 1.0e-6) {
                grid.a_y_rk4(rk_substep, j_ext) = -dt * rk41[rk_substep] * d_y_val / (k_y_val * k_y_val) / (1.0 + factor1_y);
            } else {
                grid.a_y_rk4(rk_substep, j_ext) = 0.0;
            }
            
            if (std::abs(d_y_half_val) > 1.0e-6) {
                grid.a_y_half_rk4(rk_substep, j_ext) = -dt * rk41[rk_substep] * d_y_half_val / (k_y_half_val * k_y_half_val) / (1.0 + factor1_y_half);
            } else {
                grid.a_y_half_rk4(rk_substep, j_ext) = 0.0;
            }
        }
    }
}

} // namespace seismic
