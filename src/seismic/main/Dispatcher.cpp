#include "seismic/main/Dispatcher.hpp"
#include "seismic/core/Grid.hpp"
#include "seismic/io/ModelLoader.hpp"
#include "seismic/initializers/PMLInitializer.hpp"
#include "seismic/sources/SourceTimeFunction.hpp"
#include "seismic/io/Recorder.hpp"
#include "seismic/io/Writer.hpp"
#include "seismic/solvers/TimeIntegratorFactory.hpp"
#include <iostream>
#include <cmath>

// Fortran内核声明现在位于TimeIntegrator.hpp中，这里不再需要

namespace seismic {

void run_simulation(const SimulationConfig& config) {
    std::cout << "Starting seismic simulation..." << std::endl;
    std::cout << "Grid size: " << config.grid.nx << " x " << config.grid.ny << std::endl;
    std::cout << "Time steps: " << config.time.nstep << std::endl;
    std::cout << "PML type: " << config.pml.type << std::endl;
    
    // Initialize grid with appropriate size based on PML type
    const bool use_extended_grid = (config.pml.type == "adepml");
    Grid grid(config.grid, use_extended_grid);
    
    // Load model
    std::cout << "Loading model..." << std::endl;
    load_model(grid, config.model);
    
    // Initialize PML profiles
    std::cout << "Initializing PML profiles..." << std::endl;
    initialize_pml(grid, config.pml, config.model, config.time.dt);
    
    // Create time integrator using factory pattern
    std::cout << "Creating time integrator..." << std::endl;
    try {
        auto time_integrator = create_time_integrator(config.time, config.pml);
        std::cout << "Using time integrator: " << time_integrator->get_name() << std::endl;
        std::cout << "Stability factor: " << time_integrator->get_stability_factor() << std::endl;
        
        // Check if the configuration is supported
        if (!TimeIntegratorFactory::is_supported(config.time, config.pml)) {
            std::cerr << "ERROR: Unsupported configuration!" << std::endl;
            std::cerr << TimeIntegratorFactory::get_supported_types() << std::endl;
            return;
        }
        
        // Setup seismogram recorder
        SeismogramRecorder recorder(config.receivers, config.time.nstep);
        recorder.initialize(config.grid);
        
        // Set recording mode based on PML type
        // ADE-PML uses 4-point average (original Fortran behavior)
        // CPML uses single point sampling  
        if (config.pml.type == "adepml") {
            recorder.set_recording_mode("rk4");  // 4-point average for ADE-PML
        } else {
            recorder.set_recording_mode("cpml");  // Single point for CPML
        }
        
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
        
        const Real stability_limit = time_integrator->get_stability_factor();
        
        std::cout << "Courant number: " << courant_number << std::endl;
        std::cout << "Stability limit for " << time_integrator->get_name() 
                  << ": " << stability_limit << std::endl;
                  
        if (courant_number > stability_limit) {
            std::cerr << "WARNING: Courant number > stability limit, simulation may be unstable!" << std::endl;
        }
        
        // Main time loop
        std::cout << "Starting time loop..." << std::endl;
        
        const Index display_interval = 100;
        
        for (Index it = 0; it < config.time.nstep; ++it) {
            // Execute one time step using the selected integrator
            time_integrator->step(grid, config.source, config.pml, config.time.dt, it);
            
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
                
                const Real current_time = it * config.time.dt;
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
        
    } catch (const std::runtime_error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << TimeIntegratorFactory::get_supported_types() << std::endl;
        return;
    }
}

} // namespace seismic
