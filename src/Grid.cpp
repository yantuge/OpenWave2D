#include "seismic/core/Grid.hpp"
#include <algorithm>

namespace seismic {

Grid::Grid(const GridConfig& config, bool use_extended_grid) 
    : m_nx(config.nx), m_ny(config.ny), 
      m_deltax(config.deltax), m_deltay(config.deltay),
      m_nx_extended(config.nx + 9), m_ny_extended(config.ny + 9),  // +9 for Fortran (-4:nx+4) range
      m_use_extended_grid(use_extended_grid) {
    
    const Size grid_size = static_cast<Size>(m_nx * m_ny);
    const Size extended_grid_size = static_cast<Size>(m_nx_extended * m_ny_extended);
    
    // Choose appropriate size based on PML type
    const Size actual_grid_size = use_extended_grid ? extended_grid_size : grid_size;
    const Size actual_nx = use_extended_grid ? m_nx_extended : m_nx;
    
    // Allocate wave field arrays
    m_vx.resize(actual_grid_size);
    m_vy.resize(actual_grid_size);
    m_sigmaxx.resize(actual_grid_size);
    m_sigmayy.resize(actual_grid_size);
    m_sigmaxy.resize(actual_grid_size);
    
    // Allocate material property arrays
    m_vp.resize(actual_grid_size);
    m_vs.resize(actual_grid_size);
    m_density.resize(actual_grid_size);
    m_c11.resize(actual_grid_size);
    m_c12.resize(actual_grid_size);
    m_c22.resize(actual_grid_size);
    m_c33.resize(actual_grid_size);
    m_lambda.resize(actual_grid_size);
    m_mu.resize(actual_grid_size);
    
    // Allocate RK4 intermediate arrays (3 components per field for RK4)
    const Size rk4_size = 3 * actual_grid_size;  // 3 components for RK4
    m_dvx.resize(rk4_size);
    m_dvy.resize(rk4_size);
    m_dsigmaxx.resize(rk4_size);
    m_dsigmayy.resize(rk4_size);
    m_dsigmaxy.resize(rk4_size);
    
    // Allocate PML profile arrays (1D) - conditional sizing based on PML type
    const Size profile_size_x = use_extended_grid ? m_nx_extended : m_nx;
    const Size profile_size_y = use_extended_grid ? m_ny_extended : m_ny;
    
    m_d_x.resize(profile_size_x);
    m_d_x_half.resize(profile_size_x);
    m_d_y.resize(profile_size_y);
    m_d_y_half.resize(profile_size_y);
    
    m_k_x.resize(profile_size_x);
    m_k_x_half.resize(profile_size_x);
    m_k_y.resize(profile_size_y);
    m_k_y_half.resize(profile_size_y);
    
    m_alpha_x.resize(profile_size_x);
    m_alpha_x_half.resize(profile_size_x);
    m_alpha_y.resize(profile_size_y);
    m_alpha_y_half.resize(profile_size_y);
    
    // CPML coefficients (1D)
    m_a_x.resize(profile_size_x);
    m_a_x_half.resize(profile_size_x);
    m_a_y.resize(profile_size_y);
    m_a_y_half.resize(profile_size_y);
    
    m_b_x.resize(profile_size_x);
    m_b_x_half.resize(profile_size_x);
    m_b_y.resize(profile_size_y);
    m_b_y_half.resize(profile_size_y);
    
    // ADE-PML RK4 coefficients (2D: 4 substeps × spatial points)
    m_a_x_rk4.resize(4 * m_nx_extended);
    m_a_x_half_rk4.resize(4 * m_nx_extended);
    m_a_y_rk4.resize(4 * m_ny_extended);
    m_a_y_half_rk4.resize(4 * m_ny_extended);
    
    m_b_x_rk4.resize(4 * m_nx_extended);
    m_b_x_half_rk4.resize(4 * m_nx_extended);
    m_b_y_rk4.resize(4 * m_ny_extended);
    m_b_y_half_rk4.resize(4 * m_ny_extended);
    
    // Allocate PML memory variables
    m_memory_dvx_dx.resize(actual_grid_size);
    m_memory_dvx_dy.resize(actual_grid_size);
    m_memory_dvy_dx.resize(actual_grid_size);
    m_memory_dvy_dy.resize(actual_grid_size);
    m_memory_dsigmaxx_dx.resize(actual_grid_size);
    m_memory_dsigmayy_dy.resize(actual_grid_size);
    m_memory_dsigmaxy_dx.resize(actual_grid_size);
    m_memory_dsigmaxy_dy.resize(actual_grid_size);
    
    // Initialize all fields to zero
    initialize_fields();
}

void Grid::initialize_fields() {
    // Initialize wave fields to zero
    std::fill(m_vx.begin(), m_vx.end(), 0.0);
    std::fill(m_vy.begin(), m_vy.end(), 0.0);
    std::fill(m_sigmaxx.begin(), m_sigmaxx.end(), 0.0);
    std::fill(m_sigmayy.begin(), m_sigmayy.end(), 0.0);
    std::fill(m_sigmaxy.begin(), m_sigmaxy.end(), 0.0);
    
    // Initialize PML profiles to default values
    std::fill(m_d_x.begin(), m_d_x.end(), 0.0);
    std::fill(m_d_x_half.begin(), m_d_x_half.end(), 0.0);
    std::fill(m_d_y.begin(), m_d_y.end(), 0.0);
    std::fill(m_d_y_half.begin(), m_d_y_half.end(), 0.0);
    
    std::fill(m_k_x.begin(), m_k_x.end(), 1.0);
    std::fill(m_k_x_half.begin(), m_k_x_half.end(), 1.0);
    std::fill(m_k_y.begin(), m_k_y.end(), 1.0);
    std::fill(m_k_y_half.begin(), m_k_y_half.end(), 1.0);
    
    std::fill(m_alpha_x.begin(), m_alpha_x.end(), 0.0);
    std::fill(m_alpha_x_half.begin(), m_alpha_x_half.end(), 0.0);
    std::fill(m_alpha_y.begin(), m_alpha_y.end(), 0.0);
    std::fill(m_alpha_y_half.begin(), m_alpha_y_half.end(), 0.0);
    
    std::fill(m_a_x.begin(), m_a_x.end(), 0.0);
    std::fill(m_a_x_half.begin(), m_a_x_half.end(), 0.0);
    std::fill(m_a_y.begin(), m_a_y.end(), 0.0);
    std::fill(m_a_y_half.begin(), m_a_y_half.end(), 0.0);
    
    std::fill(m_b_x.begin(), m_b_x.end(), 1.0);
    std::fill(m_b_x_half.begin(), m_b_x_half.end(), 1.0);
    std::fill(m_b_y.begin(), m_b_y.end(), 1.0);
    std::fill(m_b_y_half.begin(), m_b_y_half.end(), 1.0);
    
    // Initialize ADE-PML RK4 coefficients to default values
    std::fill(m_a_x_rk4.begin(), m_a_x_rk4.end(), 0.0);
    std::fill(m_a_x_half_rk4.begin(), m_a_x_half_rk4.end(), 0.0);
    std::fill(m_a_y_rk4.begin(), m_a_y_rk4.end(), 0.0);
    std::fill(m_a_y_half_rk4.begin(), m_a_y_half_rk4.end(), 0.0);
    
    std::fill(m_b_x_rk4.begin(), m_b_x_rk4.end(), 1.0);
    std::fill(m_b_x_half_rk4.begin(), m_b_x_half_rk4.end(), 1.0);
    std::fill(m_b_y_rk4.begin(), m_b_y_rk4.end(), 1.0);
    std::fill(m_b_y_half_rk4.begin(), m_b_y_half_rk4.end(), 1.0);
    
    // Initialize PML memory variables to zero
    std::fill(m_memory_dvx_dx.begin(), m_memory_dvx_dx.end(), 0.0);
    std::fill(m_memory_dvx_dy.begin(), m_memory_dvx_dy.end(), 0.0);
    std::fill(m_memory_dvy_dx.begin(), m_memory_dvy_dx.end(), 0.0);
    std::fill(m_memory_dvy_dy.begin(), m_memory_dvy_dy.end(), 0.0);
    std::fill(m_memory_dsigmaxx_dx.begin(), m_memory_dsigmaxx_dx.end(), 0.0);
    std::fill(m_memory_dsigmayy_dy.begin(), m_memory_dsigmayy_dy.end(), 0.0);
    std::fill(m_memory_dsigmaxy_dx.begin(), m_memory_dsigmaxy_dx.end(), 0.0);
    std::fill(m_memory_dsigmaxy_dy.begin(), m_memory_dsigmaxy_dy.end(), 0.0);
    
    // Initialize RK4 intermediate arrays to zero
    std::fill(m_dvx.begin(), m_dvx.end(), 0.0);
    std::fill(m_dvy.begin(), m_dvy.end(), 0.0);
    std::fill(m_dsigmaxx.begin(), m_dsigmaxx.end(), 0.0);
    std::fill(m_dsigmayy.begin(), m_dsigmayy.end(), 0.0);
    std::fill(m_dsigmaxy.begin(), m_dsigmaxy.end(), 0.0);
}

void Grid::backup_rk4_state() {
    // For ADE-PML RK4, we need to backup current state
    // This matches the original Fortran implementation pattern:
    // dvx(2,:,:) = vx(:,:)  and  dvx(3,:,:) = vx(:,:)
    const Size actual_grid_size = m_use_extended_grid ? 
        static_cast<Size>(m_nx_extended * m_ny_extended) :
        static_cast<Size>(m_nx * m_ny);
    
    // Backup to indices 1 and 2 (0-based: components 1 and 2) of RK4 arrays
    for (Size i = 0; i < actual_grid_size; ++i) {
        // Component 1 (dvx(2,:,:) in Fortran)
        m_dvx[actual_grid_size + i] = m_vx[i];
        m_dvy[actual_grid_size + i] = m_vy[i];
        m_dsigmaxx[actual_grid_size + i] = m_sigmaxx[i];
        m_dsigmayy[actual_grid_size + i] = m_sigmayy[i];
        m_dsigmaxy[actual_grid_size + i] = m_sigmaxy[i];
        
        // Component 2 (dvx(3,:,:) in Fortran)
        m_dvx[2 * actual_grid_size + i] = m_vx[i];
        m_dvy[2 * actual_grid_size + i] = m_vy[i];
        m_dsigmaxx[2 * actual_grid_size + i] = m_sigmaxx[i];
        m_dsigmayy[2 * actual_grid_size + i] = m_sigmayy[i];
        m_dsigmaxy[2 * actual_grid_size + i] = m_sigmaxy[i];
    }
}

} // namespace seismic
