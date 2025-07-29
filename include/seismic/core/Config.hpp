#pragma once

#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "Types.hpp"

/**
 * @file Config.hpp
 * @brief Configuration structures for seismic simulation
 */

namespace seismic {

/// Grid configuration parameters
struct GridConfig {
    Index nx;           ///< Number of grid points in X direction
    Index ny;           ///< Number of grid points in Y direction
    Real deltax;        ///< Grid spacing in X direction (meters)
    Real deltay;        ///< Grid spacing in Y direction (meters)
};

/// Physical model configuration
struct ModelConfig {
    std::string type;   ///< Model type: "homogeneous", "binary_file"
    
    // For homogeneous model
    Real vp;           ///< P-wave velocity (m/s)
    Real vs;           ///< S-wave velocity (m/s)  
    Real density;      ///< Density (kg/m³)
    
    // For anisotropic model
    Real c11, c12, c22, c33;  ///< Elastic constants
    
    // For binary file model
    std::string vp_file;      ///< P-wave velocity file path
    std::string vs_file;      ///< S-wave velocity file path
    std::string density_file; ///< Density file path
};

/// Time stepping configuration
struct TimeConfig {
    Index nstep;        ///< Total number of time steps
    Real dt;            ///< Time step size (seconds)
    std::string scheme; ///< Time scheme: "explicit", "rk4"
};

/// Source configuration
struct SourceConfig {
    Real x;             ///< Source X position (meters)
    Real y;             ///< Source Y position (meters)
    Real f0;            ///< Dominant frequency (Hz)
    Real t0;            ///< Time delay (seconds)
    Real factor;        ///< Amplitude factor
    Real angle;         ///< Force angle (degrees)
    std::string type;   ///< Source type: "ricker", "gaussian", "first_derivative"
};

/// PML absorbing boundary configuration
struct PMLConfig {
    bool use_pml_xmin;  ///< Enable PML at X minimum boundary
    bool use_pml_xmax;  ///< Enable PML at X maximum boundary
    bool use_pml_ymin;  ///< Enable PML at Y minimum boundary
    bool use_pml_ymax;  ///< Enable PML at Y maximum boundary
    Index npoints_pml;  ///< PML thickness in grid points
    Real npower;        ///< Power for PML profile
    Real k_max_pml;     ///< Maximum PML kappa value
    Real alpha_max_pml; ///< Maximum PML alpha value
    Real rcoef;         ///< Reflection coefficient
    std::string type;   ///< PML type: "cpml", "adepml"
    Real epsn;          ///< ADE-PML coefficient epsilon_n (default: 0.25)
    Real epsn1;         ///< ADE-PML coefficient epsilon_n1 (default: 0.75)
};

/// Receiver configuration
struct ReceiverConfig {
    std::vector<Real> x_positions; ///< Receiver X positions (meters)
    std::vector<Real> y_positions; ///< Receiver Y positions (meters)
};

/// Output configuration
struct OutputConfig {
    std::string output_dir;        ///< Output directory
    Index snapshot_interval;       ///< Snapshot output interval
    std::string snapshot_format;   ///< Snapshot format: "binary", "su"
    bool save_seismograms;         ///< Whether to save seismograms
    std::string seismogram_format; ///< Seismogram format: "ascii", "su"
};

/// Complete simulation configuration
struct SimulationConfig {
    GridConfig grid;
    ModelConfig model;
    TimeConfig time;
    SourceConfig source;
    PMLConfig pml;
    ReceiverConfig receivers;
    OutputConfig output;
    
    /// Load configuration from JSON
    static SimulationConfig from_json(const nlohmann::json& j);
};

} // namespace seismic
