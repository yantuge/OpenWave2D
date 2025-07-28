#include "seismic/main/Dispatcher.hpp"
#include "seismic/core/Grid.hpp"
#include "seismic/io/ModelLoader.hpp"
#include "seismic/initializers/PMLInitializer.hpp"
#include "seismic/sources/SourceTimeFunction.hpp"
#include "seismic/io/Recorder.hpp"
#include "seismic/io/Writer.hpp"
#include <iostream>
#include <cmath>

// Fortran interface declarations
extern "C" {
    void cpml_anisotropic_step_c(int nx, int ny, double deltax, double deltay, double dt, int it,
                                int isource, int jsource, double force_x, double force_y,
                                double* c11, double* c12, double* c22, double* c33, double* rho,
                                double* vx, double* vy, double* sigmaxx, double* sigmayy, double* sigmaxy,
                                double* d_x, double* d_x_half, double* d_y, double* d_y_half,
                                double* k_x, double* k_x_half, double* k_y, double* k_y_half,
                                double* alpha_x, double* alpha_x_half, double* alpha_y, double* alpha_y_half,
                                double* a_x, double* a_x_half, double* a_y, double* a_y_half,
                                double* b_x, double* b_x_half, double* b_y, double* b_y_half,
                                double* memory_dvx_dx, double* memory_dvx_dy,
                                double* memory_dvy_dx, double* memory_dvy_dy,
                                double* memory_dsigmaxx_dx, double* memory_dsigmayy_dy,
                                double* memory_dsigmaxy_dx, double* memory_dsigmaxy_dy);
    
    void adepml_rk4_step_c(int nx, int ny, double deltax, double deltay, double dt, int rk_substep,
                          int isource, int jsource, double force_x, double force_y,
                          double* lambda, double* mu, double* rho,
                          double* vx, double* vy, double* sigmaxx, double* sigmayy, double* sigmaxy,
                          double* dvx, double* dvy, double* dsigmaxx, double* dsigmayy, double* dsigmaxy,
                          double* d_x, double* d_x_half, double* d_y, double* d_y_half,
                          double* k_x, double* k_x_half, double* k_y, double* k_y_half,
                          double* alpha_x, double* alpha_x_half, double* alpha_y, double* alpha_y_half,
                          double* a_x, double* a_x_half, double* a_y, double* a_y_half,
                          double* b_x, double* b_x_half, double* b_y, double* b_y_half,
                          double* memory_dvx_dx, double* memory_dvx_dy,
                          double* memory_dvy_dx, double* memory_dvy_dy,
                          double* memory_dsigmaxx_dx, double* memory_dsigmayy_dy,
                          double* memory_dsigmaxy_dx, double* memory_dsigmaxy_dy);
}

namespace seismic {

void run_simulation(const SimulationConfig& config) {
    std::cout << "Starting seismic simulation..." << std::endl;
    std::cout << "Grid size: " << config.grid.nx << " x " << config.grid.ny << std::endl;
    std::cout << "Time steps: " << config.time.nstep << std::endl;
    std::cout << "PML type: " << config.pml.type << std::endl;
    
    // Initialize grid
    Grid grid(config.grid);
    
    // Load model
    std::cout << "Loading model..." << std::endl;
    load_model(grid, config.model);
    
    // Initialize PML profiles
    std::cout << "Initializing PML profiles..." << std::endl;
    initialize_pml(grid, config.pml, config.model, config.time.dt);
    
    // Setup seismogram recorder
    SeismogramRecorder recorder(config.receivers, config.grid, config.time.nstep);
    
    // Convert source position to grid indices
    const Index isource = static_cast<Index>(config.source.x / config.grid.deltax);
    const Index jsource = static_cast<Index>(config.source.y / config.grid.deltay);
    
    std::cout << "Source position: (" << config.source.x << ", " << config.source.y 
              << ") -> grid (" << isource << ", " << jsource << ")" << std::endl;
    
    // Check Courant stability condition
    Real cp_max = 0.0;
    if (config.model.type == "homogeneous") {
        cp_max = config.model.vp;
    } else if (config.model.type == "anisotropic") {
        cp_max = std::max(std::sqrt(config.model.c11 / config.model.density),
                         std::sqrt(config.model.c22 / config.model.density));
    } else {
        cp_max = 3000.0; // Conservative estimate
    }
    
    const Real courant_number = cp_max * config.time.dt * 
                               std::sqrt(1.0 / (config.grid.deltax * config.grid.deltax) + 
                                        1.0 / (config.grid.deltay * config.grid.deltay));
    
    std::cout << "Courant number: " << courant_number << std::endl;
    if (courant_number > 1.0) {
        std::cerr << "WARNING: Courant number > 1.0, simulation may be unstable!" << std::endl;
    }
    
    // Main time loop
    std::cout << "Starting time loop..." << std::endl;
    
    const Index display_interval = 100;
    
    for (Index it = 0; it < config.time.nstep; ++it) {
        const Real current_time = it * config.time.dt;
        
        // Calculate source amplitude and forces
        const Real source_amplitude = calculate_source_amplitude(current_time, config.source);
        Real force_x, force_y;
        calculate_source_forces(source_amplitude, config.source, force_x, force_y);
        
        // Call appropriate Fortran kernel
        if (config.pml.type == "cpml") {
            // CPML anisotropic step
            cpml_anisotropic_step_c(
                grid.nx(), grid.ny(), grid.deltax(), grid.deltay(), config.time.dt, it,
                isource, jsource, force_x, force_y,
                grid.get_c11_ptr(), grid.get_c12_ptr(), grid.get_c22_ptr(), grid.get_c33_ptr(), 
                grid.get_density_ptr(),
                grid.get_vx_ptr(), grid.get_vy_ptr(), 
                grid.get_sigmaxx_ptr(), grid.get_sigmayy_ptr(), grid.get_sigmaxy_ptr(),
                grid.get_d_x_ptr(), grid.get_d_x_half_ptr(), 
                grid.get_d_y_ptr(), grid.get_d_y_half_ptr(),
                grid.get_k_x_ptr(), grid.get_k_x_half_ptr(), 
                grid.get_k_y_ptr(), grid.get_k_y_half_ptr(),
                grid.get_alpha_x_ptr(), grid.get_alpha_x_half_ptr(), 
                grid.get_alpha_y_ptr(), grid.get_alpha_y_half_ptr(),
                grid.get_a_x_ptr(), grid.get_a_x_half_ptr(), 
                grid.get_a_y_ptr(), grid.get_a_y_half_ptr(),
                grid.get_b_x_ptr(), grid.get_b_x_half_ptr(), 
                grid.get_b_y_ptr(), grid.get_b_y_half_ptr(),
                grid.get_memory_dvx_dx_ptr(), grid.get_memory_dvx_dy_ptr(),
                grid.get_memory_dvy_dx_ptr(), grid.get_memory_dvy_dy_ptr(),
                grid.get_memory_dsigmaxx_dx_ptr(), grid.get_memory_dsigmayy_dy_ptr(),
                grid.get_memory_dsigmaxy_dx_ptr(), grid.get_memory_dsigmaxy_dy_ptr()
            );
        }
        else if (config.pml.type == "adepml") {
            // ADE-PML RK4 implementation - backup current state for RK4
            grid.backup_rk4_state();
            
            // RK4 loop (4 substeps) - following original Fortran implementation
            for (int rk_substep = 1; rk_substep <= 4; ++rk_substep) {
                // Calculate source forces for this RK4 substep
                const Real rk41_coeffs[4] = {0.0, 0.5, 0.5, 1.0};
                const Real substep_time = current_time + rk41_coeffs[rk_substep-1] * config.time.dt;
                const Real substep_source_amplitude = calculate_source_amplitude(substep_time, config.source);
                Real substep_force_x, substep_force_y;
                calculate_source_forces(substep_source_amplitude, config.source, substep_force_x, substep_force_y);
                
                // Call ADE-PML RK4 Fortran kernel
                adepml_rk4_step_c(
                    grid.nx(), grid.ny(), grid.deltax(), grid.deltay(), config.time.dt, rk_substep,
                    isource, jsource, substep_force_x, substep_force_y,
                    grid.get_lambda_ptr(), grid.get_mu_ptr(), grid.get_density_ptr(),
                    grid.get_vx_ptr(), grid.get_vy_ptr(), 
                    grid.get_sigmaxx_ptr(), grid.get_sigmayy_ptr(), grid.get_sigmaxy_ptr(),
                    grid.get_dvx_ptr(), grid.get_dvy_ptr(), 
                    grid.get_dsigmaxx_ptr(), grid.get_dsigmayy_ptr(), grid.get_dsigmaxy_ptr(),
                    grid.get_d_x_ptr(), grid.get_d_x_half_ptr(), 
                    grid.get_d_y_ptr(), grid.get_d_y_half_ptr(),
                    grid.get_k_x_ptr(), grid.get_k_x_half_ptr(), 
                    grid.get_k_y_ptr(), grid.get_k_y_half_ptr(),
                    grid.get_alpha_x_ptr(), grid.get_alpha_x_half_ptr(), 
                    grid.get_alpha_y_ptr(), grid.get_alpha_y_half_ptr(),
                    grid.get_a_x_ptr(), grid.get_a_x_half_ptr(), 
                    grid.get_a_y_ptr(), grid.get_a_y_half_ptr(),
                    grid.get_b_x_ptr(), grid.get_b_x_half_ptr(), 
                    grid.get_b_y_ptr(), grid.get_b_y_half_ptr(),
                    grid.get_memory_dvx_dx_ptr(), grid.get_memory_dvx_dy_ptr(),
                    grid.get_memory_dvy_dx_ptr(), grid.get_memory_dvy_dy_ptr(),
                    grid.get_memory_dsigmaxx_dx_ptr(), grid.get_memory_dsigmayy_dy_ptr(),
                    grid.get_memory_dsigmaxy_dx_ptr(), grid.get_memory_dsigmaxy_dy_ptr()
                );
            }
            // End of RK4 loop
        }
        else {
            std::cerr << "ERROR: Unknown PML type: " << config.pml.type << std::endl;
            return;
        }
        
        // Record seismograms
        recorder.record_time_step(grid, it);
        
        // Save snapshots
        if (config.output.snapshot_interval > 0 && it % config.output.snapshot_interval == 0) {
            save_snapshot(grid, config.output, it, "vx");
            save_snapshot(grid, config.output, it, "vy");
        }
        
        // Display progress
        if (it % display_interval == 0 || it == config.time.nstep - 1) {
            // Calculate maximum velocity for stability check
            Real max_velocity = 0.0;
            for (Index j = 0; j < grid.ny(); ++j) {
                for (Index i = 0; i < grid.nx(); ++i) {
                    const Real vel_magnitude = std::sqrt(grid.vx(i,j) * grid.vx(i,j) + 
                                                         grid.vy(i,j) * grid.vy(i,j));
                    max_velocity = std::max(max_velocity, vel_magnitude);
                }
            }
            
            std::cout << "Time step " << it << "/" << config.time.nstep 
                      << ", Time: " << current_time << " s"
                      << ", Max velocity: " << max_velocity << " m/s" << std::endl;
            
            // Check for instability
            if (max_velocity > 1.0e25) {
                std::cerr << "ERROR: Simulation became unstable!" << std::endl;
                break;
            }
        }
    }
    
    // Save final results
    std::cout << "Saving results..." << std::endl;
    recorder.save_seismograms(config.output, config.time.dt);
    
    std::cout << "Simulation completed successfully!" << std::endl;
}

} // namespace seismic
