#include "seismic/io/ModelLoader.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace seismic {

void load_model(Grid& grid, const ModelConfig& config) {
    if (config.type == "homogeneous") {
        build_homogeneous_model(grid, config);
    } else if (config.type == "anisotropic") {
        build_anisotropic_model(grid, config);
    } else if (config.type == "binary_file") {
        load_model_from_binary(grid, config);
    } else {
        throw std::runtime_error("Unknown model type: " + config.type);
    }
}

void build_homogeneous_model(Grid& grid, const ModelConfig& config) {
    const Size grid_size = static_cast<Size>(grid.nx() * grid.ny());
    
    // Fill velocity and density arrays with constant values
    std::fill_n(grid.get_vp_ptr(), grid_size, config.vp);
    std::fill_n(grid.get_vs_ptr(), grid_size, config.vs);
    std::fill_n(grid.get_density_ptr(), grid_size, config.density);
    
    // Compute elastic constants for isotropic case
    const Real lambda = config.density * (config.vp * config.vp - 2.0 * config.vs * config.vs);
    const Real mu = config.density * config.vs * config.vs;
    
    // Fill elastic constant arrays
    std::fill_n(grid.get_c11_ptr(), grid_size, lambda + 2.0 * mu);  // c11
    std::fill_n(grid.get_c12_ptr(), grid_size, lambda);             // c12
    std::fill_n(grid.get_c22_ptr(), grid_size, lambda + 2.0 * mu);  // c22 
    std::fill_n(grid.get_c33_ptr(), grid_size, mu);                 // c33
}

void build_anisotropic_model(Grid& grid, const ModelConfig& config) {
    const Size grid_size = static_cast<Size>(grid.nx() * grid.ny());
    
    // Fill with anisotropic elastic constants
    std::fill_n(grid.get_c11_ptr(), grid_size, config.c11);
    std::fill_n(grid.get_c12_ptr(), grid_size, config.c12);
    std::fill_n(grid.get_c22_ptr(), grid_size, config.c22);
    std::fill_n(grid.get_c33_ptr(), grid_size, config.c33);
    std::fill_n(grid.get_density_ptr(), grid_size, config.density);
    
    // Compute equivalent velocities for PML initialization
    const Real vp_horizontal = std::sqrt(config.c11 / config.density);
    const Real vp_vertical = std::sqrt(config.c22 / config.density);
    const Real vs = std::sqrt(config.c33 / config.density);
    
    std::fill_n(grid.get_vp_ptr(), grid_size, std::max(vp_horizontal, vp_vertical));
    std::fill_n(grid.get_vs_ptr(), grid_size, vs);
}

void load_model_from_binary(Grid& grid, const ModelConfig& config) {
    const Size grid_size = static_cast<Size>(grid.nx() * grid.ny());
    
    // Load P-wave velocity
    if (!config.vp_file.empty()) {
        std::ifstream vp_stream(config.vp_file, std::ios::binary);
        if (!vp_stream) {
            throw std::runtime_error("Cannot open P-wave velocity file: " + config.vp_file);
        }
        vp_stream.read(reinterpret_cast<char*>(grid.get_vp_ptr()), 
                      grid_size * sizeof(Real));
        vp_stream.close();
    }
    
    // Load S-wave velocity
    if (!config.vs_file.empty()) {
        std::ifstream vs_stream(config.vs_file, std::ios::binary);
        if (!vs_stream) {
            throw std::runtime_error("Cannot open S-wave velocity file: " + config.vs_file);
        }
        vs_stream.read(reinterpret_cast<char*>(grid.get_vs_ptr()), 
                      grid_size * sizeof(Real));
        vs_stream.close();
    }
    
    // Load density
    if (!config.density_file.empty()) {
        std::ifstream density_stream(config.density_file, std::ios::binary);
        if (!density_stream) {
            throw std::runtime_error("Cannot open density file: " + config.density_file);
        }
        density_stream.read(reinterpret_cast<char*>(grid.get_density_ptr()), 
                           grid_size * sizeof(Real));
        density_stream.close();
    }
    
    // Compute elastic constants from velocities and density
    for (Index j = 0; j < grid.ny(); ++j) {
        for (Index i = 0; i < grid.nx(); ++i) {
            const Real rho = grid.density(i, j);
            const Real vp = grid.vp(i, j);
            const Real vs = grid.vs(i, j);
            
            const Real lambda = rho * (vp * vp - 2.0 * vs * vs);
            const Real mu = rho * vs * vs;
            
            grid.c11(i, j) = lambda + 2.0 * mu;
            grid.c12(i, j) = lambda;
            grid.c22(i, j) = lambda + 2.0 * mu;
            grid.c33(i, j) = mu;
        }
    }
}

} // namespace seismic
