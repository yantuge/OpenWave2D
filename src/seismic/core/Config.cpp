#include "seismic/core/Config.hpp"
#include <stdexcept>
#include <cmath>

namespace seismic {

SimulationConfig SimulationConfig::from_json(const nlohmann::json& j) {
    SimulationConfig config;
    
    // Parse grid configuration
    if (j.contains("grid")) {
        const auto& grid_json = j["grid"];
        config.grid.nx = grid_json.value("nx", 401);
        config.grid.ny = grid_json.value("ny", 401);
        config.grid.deltax = grid_json.value("deltax", 0.0625e-2);
        config.grid.deltay = grid_json.value("deltay", config.grid.deltax);
    }
    
    // Parse model configuration
    if (j.contains("model")) {
        const auto& model_json = j["model"];
        config.model.type = model_json.value("type", "homogeneous");
        
        if (config.model.type == "homogeneous") {
            config.model.vp = model_json.value("vp", 2000.0);
            config.model.vs = model_json.value("vs", 1150.0);
            config.model.density = model_json.value("density", 2000.0);
        } else if (config.model.type == "anisotropic") {
            config.model.c11 = model_json.value("c11", 4.0e10);
            config.model.c12 = model_json.value("c12", 3.8e10);
            config.model.c22 = model_json.value("c22", 20.0e10);
            config.model.c33 = model_json.value("c33", 2.0e10);
            config.model.density = model_json.value("density", 4000.0);
        } else if (config.model.type == "binary_file") {
            config.model.vp_file = model_json.value("vp_file", "");
            config.model.vs_file = model_json.value("vs_file", "");
            config.model.density_file = model_json.value("density_file", "");
        }
    }
    
    // Parse time configuration
    if (j.contains("time")) {
        const auto& time_json = j["time"];
        config.time.nstep = time_json.value("nstep", 3000);
        config.time.dt = time_json.value("dt", 50.0e-9);
        config.time.scheme = time_json.value("scheme", "explicit");
    }
    
    // Parse source configuration
    if (j.contains("source")) {
        const auto& source_json = j["source"];
        config.source.x = source_json.value("x", config.grid.nx * config.grid.deltax / 2.0);
        config.source.y = source_json.value("y", config.grid.ny * config.grid.deltay / 2.0);
        config.source.f0 = source_json.value("f0", 200000.0);
        config.source.t0 = source_json.value("t0", 1.2 / config.source.f0);
        config.source.factor = source_json.value("factor", 1.0e7);
        config.source.angle = source_json.value("angle", 0.0);
        config.source.type = source_json.value("type", "ricker");
    }
    
    // Parse PML configuration
    if (j.contains("pml")) {
        const auto& pml_json = j["pml"];
        config.pml.use_pml_xmin = pml_json.value("use_pml_xmin", true);
        config.pml.use_pml_xmax = pml_json.value("use_pml_xmax", true);
        config.pml.use_pml_ymin = pml_json.value("use_pml_ymin", true);
        config.pml.use_pml_ymax = pml_json.value("use_pml_ymax", true);
        config.pml.npoints_pml = pml_json.value("npoints_pml", 10);
        config.pml.npower = pml_json.value("npower", 2.0);
        config.pml.k_max_pml = pml_json.value("k_max_pml", 1.0);
        config.pml.alpha_max_pml = pml_json.value("alpha_max_pml", 
                                                  2.0 * 3.141592653589793 * (config.source.f0 / 2.0));
        config.pml.rcoef = pml_json.value("rcoef", 0.001);
        config.pml.type = pml_json.value("type", "cpml");
        // ADE-PML specific parameters
        config.pml.epsn = pml_json.value("epsn", 0.25);
        config.pml.epsn1 = pml_json.value("epsn1", 0.75);
    }
    
    // Parse receiver configuration
    if (j.contains("receivers")) {
        const auto& receivers_json = j["receivers"];
        if (receivers_json.is_array()) {
            for (const auto& receiver : receivers_json) {
                config.receivers.x_positions.push_back(receiver["x"]);
                config.receivers.y_positions.push_back(receiver["y"]);
            }
        }
    }
    
    // Parse output configuration
    if (j.contains("output")) {
        const auto& output_json = j["output"];
        config.output.output_dir = output_json.value("output_dir", "./output");
        config.output.snapshot_interval = output_json.value("snapshot_interval", 100);
        config.output.snapshot_format = output_json.value("snapshot_format", "binary");
        config.output.save_seismograms = output_json.value("save_seismograms", true);
        config.output.seismogram_format = output_json.value("seismogram_format", "su");
    }
    
    return config;
}

} // namespace seismic
