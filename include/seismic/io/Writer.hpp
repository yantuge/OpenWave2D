#pragma once

#include "seismic/core/Types.hpp"
#include "seismic/core/Config.hpp"
#include "seismic/core/Grid.hpp"

/**
 * @file Writer.hpp
 * @brief Snapshot and data writing functionality
 */

namespace seismic {

/**
 * @brief Save velocity snapshot to binary file
 * @param grid Current grid state
 * @param output_config Output configuration
 * @param it Current time step
 * @param component Component to save ("vx" or "vy")
 */
void save_snapshot(const Grid& grid, const OutputConfig& output_config, 
                   Index it, const std::string& component);

/**
 * @brief Save velocity snapshot to SU format
 * @param grid Current grid state
 * @param output_config Output configuration
 * @param it Current time step
 * @param component Component to save ("vx" or "vy")
 */
void save_snapshot_su(const Grid& grid, const OutputConfig& output_config,
                      Index it, const std::string& component);

/**
 * @brief Save velocity snapshot to SU format matching original Fortran code
 * @param grid Current grid state
 * @param output_config Output configuration
 * @param it Current time step
 * @param component Component to save ("vx" or "vy")
 */
void save_snapshot_su_fortran_style(const Grid& grid, const OutputConfig& output_config,
                                    Index it, const std::string& component);

/**
 * @brief Write seismogram data in ASCII format
 * @param seismograms 2D vector of seismogram data [receiver][time]
 * @param filename Output filename
 * @param dt Time step size
 */
void write_seismograms_ascii(const std::vector<std::vector<Real>>& seismograms,
                            const std::string& filename, Real dt);

/**
 * @brief Write seismogram data in SU format
 * @param seismograms 2D vector of seismogram data [receiver][time]
 * @param filename Output filename
 * @param dt Time step size
 */
void write_seismograms_su(const std::vector<std::vector<Real>>& seismograms,
                         const std::string& filename, Real dt);

} // namespace seismic
