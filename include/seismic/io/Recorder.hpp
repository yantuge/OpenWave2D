#pragma once

#include <vector>
#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"
#include "seismic/core/Grid.hpp"

/**
 * @file Recorder.hpp
 * @brief Seismogram recording functionality
 */

namespace seismic {

/**
 * @class SeismogramRecorder
 * @brief Records seismograms at receiver locations
 */
class SeismogramRecorder {
public:
    /**
     * @brief Constructor
     * @param receivers Receiver configuration
     * @param grid_config Grid configuration for coordinate conversion
     * @param nstep Total number of time steps
     */
    SeismogramRecorder(const ReceiverConfig& receivers, 
                      const GridConfig& grid_config,
                      Index nstep);

    /**
     * @brief Record current time step
     * @param grid Current grid state
     * @param it Current time step index
     */
    void record_time_step(const Grid& grid, Index it);

    /**
     * @brief Save seismograms to files
     * @param output_config Output configuration
     * @param dt Time step size
     */
    void save_seismograms(const OutputConfig& output_config, Real dt) const;

private:
    struct ReceiverLocation {
        Index i, j;    ///< Grid indices
        Real x, y;     ///< Physical coordinates
    };

    std::vector<ReceiverLocation> m_receivers;
    std::vector<std::vector<Real>> m_seismograms_vx;  ///< Vx seismograms
    std::vector<std::vector<Real>> m_seismograms_vy;  ///< Vy seismograms
    Index m_nstep;
    
    /// Convert physical coordinates to grid indices
    void compute_receiver_locations(const ReceiverConfig& receivers, 
                                   const GridConfig& grid_config);
};

} // namespace seismic
