#include "seismic/io/Recorder.hpp"
#include "seismic/io/Writer.hpp"
#include <cmath>
#include <algorithm>

namespace seismic {

SeismogramRecorder::SeismogramRecorder(const ReceiverConfig& receivers, 
                                      const GridConfig& grid_config,
                                      Index nstep) 
    : m_nstep(nstep) {
    compute_receiver_locations(receivers, grid_config);
    
    const Size num_receivers = m_receivers.size();
    m_seismograms_vx.resize(num_receivers);
    m_seismograms_vy.resize(num_receivers);
    
    for (Size i = 0; i < num_receivers; ++i) {
        m_seismograms_vx[i].resize(nstep, 0.0);
        m_seismograms_vy[i].resize(nstep, 0.0);
    }
}

void SeismogramRecorder::compute_receiver_locations(const ReceiverConfig& receivers, 
                                                   const GridConfig& grid_config) {
    const Size num_receivers = receivers.x_positions.size();
    m_receivers.resize(num_receivers);
    
    for (Size r = 0; r < num_receivers; ++r) {
        const Real x_target = receivers.x_positions[r];
        const Real y_target = receivers.y_positions[r];
        
        // Find closest grid point
        Real min_distance = 1.0e30;
        Index best_i = 0, best_j = 0;
        
        for (Index j = 0; j < grid_config.ny; ++j) {
            for (Index i = 0; i < grid_config.nx; ++i) {
                const Real x_grid = i * grid_config.deltax;
                const Real y_grid = j * grid_config.deltay;
                
                const Real distance = std::sqrt((x_grid - x_target) * (x_grid - x_target) +
                                               (y_grid - y_target) * (y_grid - y_target));
                
                if (distance < min_distance) {
                    min_distance = distance;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        
        m_receivers[r].i = best_i;
        m_receivers[r].j = best_j;
        m_receivers[r].x = best_i * grid_config.deltax;
        m_receivers[r].y = best_j * grid_config.deltay;
    }
}

void SeismogramRecorder::record_time_step(const Grid& grid, Index it) {
    for (Size r = 0; r < m_receivers.size(); ++r) {
        const Index i = m_receivers[r].i;
        const Index j = m_receivers[r].j;
        
        if (it < m_nstep) {
            m_seismograms_vx[r][it] = grid.vx(i, j);
            m_seismograms_vy[r][it] = grid.vy(i, j);
        }
    }
}

void SeismogramRecorder::save_seismograms(const OutputConfig& output_config, Real dt) const {
    if (!output_config.save_seismograms) return;
    
    const std::string vx_filename = output_config.output_dir + "/seismograms_vx." + output_config.seismogram_format;
    const std::string vy_filename = output_config.output_dir + "/seismograms_vy." + output_config.seismogram_format;
    
    if (output_config.seismogram_format == "ascii") {
        write_seismograms_ascii(m_seismograms_vx, vx_filename, dt);
        write_seismograms_ascii(m_seismograms_vy, vy_filename, dt);
    } else if (output_config.seismogram_format == "su") {
        write_seismograms_su(m_seismograms_vx, vx_filename, dt);
        write_seismograms_su(m_seismograms_vy, vy_filename, dt);
    }
}

} // namespace seismic
