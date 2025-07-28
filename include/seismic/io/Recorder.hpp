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
     * @param receiver_config Single receiver configuration
     * @param nstep Total number of time steps
     */
    SeismogramRecorder(const ReceiverConfig& receiver_config, Index nstep);

    /**
     * @brief Initialize receiver locations based on grid configuration
     * @param grid_config Grid configuration for coordinate conversion
     */
    void initialize(const GridConfig& grid_config);

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
    /**
     * @brief Save seismograms in Fortran-style ASCII format (one file per receiver)
     * @param output_config Output configuration
     * @param dt Time step size
     */
    void save_seismograms_fortran_style(const OutputConfig& output_config, Real dt) const;

    /**
     * @brief Write seismograms in ASCII format
     * @param seismograms Seismogram data
     * @param filename Output filename
     * @param dt Time step size
     */
    void write_seismograms_ascii(const std::vector<std::vector<Real>>& seismograms, 
                                const std::string& filename, Real dt) const;

    /**
     * @brief Write seismograms in SU binary format
     * @param seismograms Seismogram data
     * @param filename Output filename
     * @param dt Time step size
     */
    void write_seismograms_su(const std::vector<std::vector<Real>>& seismograms, 
                             const std::string& filename, Real dt) const;

private:
    struct ReceiverLocation {
        Index i, j;    ///< Grid indices
        Real x, y;     ///< Physical coordinates
    };

    std::vector<ReceiverLocation> m_receivers;
    std::vector<std::vector<Real>> m_seismograms_vx;  ///< Vx seismograms
    std::vector<std::vector<Real>> m_seismograms_vy;  ///< Vy seismograms
    Index m_nstep;
};

} // namespace seismic
