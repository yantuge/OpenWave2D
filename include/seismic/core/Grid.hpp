#pragma once

#include <vector>
#include <memory>
#include "Types.hpp"
#include "Config.hpp"

/**
 * @file Grid.hpp
 * @brief Grid class for managing 2D seismic simulation data
 */

namespace seismic {

/**
 * @class Grid
 * @brief Manages all 2D field arrays for seismic simulation
 * 
 * The Grid class uses column-major ordering for efficient Fortran interoperability.
 * All arrays are stored as 1D vectors and accessed through 2D accessors.
 */
class Grid {
public:
    /**
     * @brief Construct grid from configuration
     * @param config Grid configuration parameters
     * @param use_extended_grid Whether to use extended grid for ADE-PML
     */
    explicit Grid(const GridConfig& config, bool use_extended_grid = false);

    /// Destructor
    ~Grid() = default;

    // Grid dimensions accessors
    Index nx() const { return m_nx; }
    Index ny() const { return m_ny; }
    Real deltax() const { return m_deltax; }
    Real deltay() const { return m_deltay; }
    
    // For Fortran kernels: get the grid dimensions that should be passed to Fortran
    // For ADE-PML: nx (because Fortran expects (-4:nx+4) and we have nx+9 total)
    // For CPML: nx (regular grid)
    Index fortran_nx() const { 
        return m_nx; 
    }
    Index fortran_ny() const { 
        return m_ny; 
    }

    // Velocity field accessors (automatically choose correct grid size)
    Real& vx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;  // Convert from regular grid to extended grid
            const Index j_ext = j + 4;
            return m_vx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_vx[j * m_nx + i]; 
        }
    }
    const Real& vx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_vx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_vx[j * m_nx + i]; 
        }
    }
    
    Real& vy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_vy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_vy[j * m_nx + i]; 
        }
    }
    const Real& vy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_vy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_vy[j * m_nx + i]; 
        }
    }
    
    // Velocity field accessors for extended grid (ADE-PML with ghost points)
    Real& vx_extended(Index i, Index j) { 
        const Index i_ext = i + 4;  // Convert from regular grid to extended grid
        const Index j_ext = j + 4;
        return m_vx[j_ext * m_nx_extended + i_ext]; 
    }
    const Real& vx_extended(Index i, Index j) const { 
        const Index i_ext = i + 4;
        const Index j_ext = j + 4;
        return m_vx[j_ext * m_nx_extended + i_ext]; 
    }
    
    Real& vy_extended(Index i, Index j) { 
        const Index i_ext = i + 4;
        const Index j_ext = j + 4;
        return m_vy[j_ext * m_nx_extended + i_ext]; 
    }
    const Real& vy_extended(Index i, Index j) const { 
        const Index i_ext = i + 4;
        const Index j_ext = j + 4;
        return m_vy[j_ext * m_nx_extended + i_ext]; 
    }

    // Stress field accessors (automatically choose correct grid size)
    Real& sigmaxx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmaxx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmaxx[j * m_nx + i]; 
        }
    }
    const Real& sigmaxx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmaxx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmaxx[j * m_nx + i]; 
        }
    }
    
    Real& sigmayy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmayy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmayy[j * m_nx + i]; 
        }
    }
    const Real& sigmayy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmayy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmayy[j * m_nx + i]; 
        }
    }
    
    Real& sigmaxy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmaxy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmaxy[j * m_nx + i]; 
        }
    }
    const Real& sigmaxy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_sigmaxy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_sigmaxy[j * m_nx + i]; 
        }
    }

    // Material property accessors
    Real& vp(Index i, Index j) { return m_vp[j * m_nx + i]; }
    const Real& vp(Index i, Index j) const { return m_vp[j * m_nx + i]; }
    
    Real& vs(Index i, Index j) { return m_vs[j * m_nx + i]; }
    const Real& vs(Index i, Index j) const { return m_vs[j * m_nx + i]; }
    
    Real& density(Index i, Index j) { return m_density[j * m_nx + i]; }
    const Real& density(Index i, Index j) const { return m_density[j * m_nx + i]; }

    // Anisotropic material properties
    Real& c11(Index i, Index j) { return m_c11[j * m_nx + i]; }
    const Real& c11(Index i, Index j) const { return m_c11[j * m_nx + i]; }
    
    Real& c12(Index i, Index j) { return m_c12[j * m_nx + i]; }
    const Real& c12(Index i, Index j) const { return m_c12[j * m_nx + i]; }
    
    Real& c22(Index i, Index j) { return m_c22[j * m_nx + i]; }
    const Real& c22(Index i, Index j) const { return m_c22[j * m_nx + i]; }
    
    Real& c33(Index i, Index j) { return m_c33[j * m_nx + i]; }
    const Real& c33(Index i, Index j) const { return m_c33[j * m_nx + i]; }

    // PML profile accessors (1D arrays)
    Real& d_x(Index i) { return m_d_x[i]; }
    const Real& d_x(Index i) const { return m_d_x[i]; }
    
    Real& d_x_half(Index i) { return m_d_x_half[i]; }
    const Real& d_x_half(Index i) const { return m_d_x_half[i]; }
    
    Real& d_y(Index j) { return m_d_y[j]; }
    const Real& d_y(Index j) const { return m_d_y[j]; }
    
    Real& d_y_half(Index j) { return m_d_y_half[j]; }
    const Real& d_y_half(Index j) const { return m_d_y_half[j]; }

    Real& k_x(Index i) { return m_k_x[i]; }
    const Real& k_x(Index i) const { return m_k_x[i]; }
    
    Real& k_x_half(Index i) { return m_k_x_half[i]; }
    const Real& k_x_half(Index i) const { return m_k_x_half[i]; }
    
    Real& k_y(Index j) { return m_k_y[j]; }
    const Real& k_y(Index j) const { return m_k_y[j]; }
    
    Real& k_y_half(Index j) { return m_k_y_half[j]; }
    const Real& k_y_half(Index j) const { return m_k_y_half[j]; }

    Real& alpha_x(Index i) { return m_alpha_x[i]; }
    const Real& alpha_x(Index i) const { return m_alpha_x[i]; }
    
    Real& alpha_x_half(Index i) { return m_alpha_x_half[i]; }
    const Real& alpha_x_half(Index i) const { return m_alpha_x_half[i]; }
    
    Real& alpha_y(Index j) { return m_alpha_y[j]; }
    const Real& alpha_y(Index j) const { return m_alpha_y[j]; }
    
    Real& alpha_y_half(Index j) { return m_alpha_y_half[j]; }
    const Real& alpha_y_half(Index j) const { return m_alpha_y_half[j]; }

    // For CPML (simple 1D access)
    Real& a_x(Index i) { return m_a_x[i]; }
    const Real& a_x(Index i) const { return m_a_x[i]; }
    
    Real& a_x_half(Index i) { return m_a_x_half[i]; }
    const Real& a_x_half(Index i) const { return m_a_x_half[i]; }
    
    Real& a_y(Index j) { return m_a_y[j]; }
    const Real& a_y(Index j) const { return m_a_y[j]; }
    
    Real& a_y_half(Index j) { return m_a_y_half[j]; }
    const Real& a_y_half(Index j) const { return m_a_y_half[j]; }

    Real& b_x(Index i) { return m_b_x[i]; }
    const Real& b_x(Index i) const { return m_b_x[i]; }
    
    Real& b_x_half(Index i) { return m_b_x_half[i]; }
    const Real& b_x_half(Index i) const { return m_b_x_half[i]; }
    
    Real& b_y(Index j) { return m_b_y[j]; }
    const Real& b_y(Index j) const { return m_b_y[j]; }
    
    Real& b_y_half(Index j) { return m_b_y_half[j]; }
    const Real& b_y_half(Index j) const { return m_b_y_half[j]; }

    // For ADE-PML (2D access: [rk_substep][spatial_index])
    Real& a_x_rk4(Index rk_substep, Index i) { return m_a_x_rk4[rk_substep * m_nx_extended + i]; }
    const Real& a_x_rk4(Index rk_substep, Index i) const { return m_a_x_rk4[rk_substep * m_nx_extended + i]; }
    
    Real& a_x_half_rk4(Index rk_substep, Index i) { return m_a_x_half_rk4[rk_substep * m_nx_extended + i]; }
    const Real& a_x_half_rk4(Index rk_substep, Index i) const { return m_a_x_half_rk4[rk_substep * m_nx_extended + i]; }
    
    Real& a_y_rk4(Index rk_substep, Index j) { return m_a_y_rk4[rk_substep * m_ny_extended + j]; }
    const Real& a_y_rk4(Index rk_substep, Index j) const { return m_a_y_rk4[rk_substep * m_ny_extended + j]; }
    
    Real& a_y_half_rk4(Index rk_substep, Index j) { return m_a_y_half_rk4[rk_substep * m_ny_extended + j]; }
    const Real& a_y_half_rk4(Index rk_substep, Index j) const { return m_a_y_half_rk4[rk_substep * m_ny_extended + j]; }

    Real& b_x_rk4(Index rk_substep, Index i) { return m_b_x_rk4[rk_substep * m_nx_extended + i]; }
    const Real& b_x_rk4(Index rk_substep, Index i) const { return m_b_x_rk4[rk_substep * m_nx_extended + i]; }
    
    Real& b_x_half_rk4(Index rk_substep, Index i) { return m_b_x_half_rk4[rk_substep * m_nx_extended + i]; }
    const Real& b_x_half_rk4(Index rk_substep, Index i) const { return m_b_x_half_rk4[rk_substep * m_nx_extended + i]; }
    
    Real& b_y_rk4(Index rk_substep, Index j) { return m_b_y_rk4[rk_substep * m_ny_extended + j]; }
    const Real& b_y_rk4(Index rk_substep, Index j) const { return m_b_y_rk4[rk_substep * m_ny_extended + j]; }
    
    Real& b_y_half_rk4(Index rk_substep, Index j) { return m_b_y_half_rk4[rk_substep * m_ny_extended + j]; }
    const Real& b_y_half_rk4(Index rk_substep, Index j) const { return m_b_y_half_rk4[rk_substep * m_ny_extended + j]; }

    // Memory variables for PML (2D arrays) - automatically choose correct grid size
    Real& memory_dvx_dx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvx_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvx_dx[j * m_nx + i]; 
        }
    }
    const Real& memory_dvx_dx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvx_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvx_dx[j * m_nx + i]; 
        }
    }
    
    Real& memory_dvx_dy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvx_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvx_dy[j * m_nx + i]; 
        }
    }
    const Real& memory_dvx_dy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvx_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvx_dy[j * m_nx + i]; 
        }
    }
    
    Real& memory_dvy_dx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvy_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvy_dx[j * m_nx + i]; 
        }
    }
    const Real& memory_dvy_dx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvy_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvy_dx[j * m_nx + i]; 
        }
    }
    
    Real& memory_dvy_dy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvy_dy[j * m_nx + i]; 
        }
    }
    const Real& memory_dvy_dy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dvy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dvy_dy[j * m_nx + i]; 
        }
    }
    
    Real& memory_dsigmaxx_dx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxx_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxx_dx[j * m_nx + i]; 
        }
    }
    const Real& memory_dsigmaxx_dx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxx_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxx_dx[j * m_nx + i]; 
        }
    }
    
    Real& memory_dsigmayy_dy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmayy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmayy_dy[j * m_nx + i]; 
        }
    }
    const Real& memory_dsigmayy_dy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmayy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmayy_dy[j * m_nx + i]; 
        }
    }
    
    Real& memory_dsigmaxy_dx(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxy_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxy_dx[j * m_nx + i]; 
        }
    }
    const Real& memory_dsigmaxy_dx(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxy_dx[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxy_dx[j * m_nx + i]; 
        }
    }
    
    Real& memory_dsigmaxy_dy(Index i, Index j) { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxy_dy[j * m_nx + i]; 
        }
    }
    const Real& memory_dsigmaxy_dy(Index i, Index j) const { 
        if (m_use_extended_grid) {
            const Index i_ext = i + 4;
            const Index j_ext = j + 4;
            return m_memory_dsigmaxy_dy[j_ext * m_nx_extended + i_ext]; 
        } else {
            return m_memory_dsigmaxy_dy[j * m_nx + i]; 
        }
    }

    // Raw pointer accessors for Fortran interface
    Real* get_vx_ptr() { return m_vx.data(); }
    Real* get_vy_ptr() { return m_vy.data(); }
    const Real* get_vx_ptr() const { return m_vx.data(); }
    const Real* get_vy_ptr() const { return m_vy.data(); }
    Real* get_sigmaxx_ptr() { return m_sigmaxx.data(); }
    Real* get_sigmayy_ptr() { return m_sigmayy.data(); }
    Real* get_sigmaxy_ptr() { return m_sigmaxy.data(); }
    
    Real* get_vp_ptr() { return m_vp.data(); }
    Real* get_vs_ptr() { return m_vs.data(); }
    Real* get_density_ptr() { return m_density.data(); }
    
    Real* get_c11_ptr() { return m_c11.data(); }
    Real* get_c12_ptr() { return m_c12.data(); }
    Real* get_c22_ptr() { return m_c22.data(); }
    Real* get_c33_ptr() { return m_c33.data(); }
    
    Real* get_d_x_ptr() { return m_d_x.data(); }
    Real* get_d_x_half_ptr() { return m_d_x_half.data(); }
    Real* get_d_y_ptr() { return m_d_y.data(); }
    Real* get_d_y_half_ptr() { return m_d_y_half.data(); }
    
    Real* get_k_x_ptr() { return m_k_x.data(); }
    Real* get_k_x_half_ptr() { return m_k_x_half.data(); }
    Real* get_k_y_ptr() { return m_k_y.data(); }
    Real* get_k_y_half_ptr() { return m_k_y_half.data(); }
    
    Real* get_alpha_x_ptr() { return m_alpha_x.data(); }
    Real* get_alpha_x_half_ptr() { return m_alpha_x_half.data(); }
    Real* get_alpha_y_ptr() { return m_alpha_y.data(); }
    Real* get_alpha_y_half_ptr() { return m_alpha_y_half.data(); }
    
    Real* get_a_x_ptr() { return m_a_x.data(); }
    Real* get_a_x_half_ptr() { return m_a_x_half.data(); }
    Real* get_a_y_ptr() { return m_a_y.data(); }
    Real* get_a_y_half_ptr() { return m_a_y_half.data(); }
    
    Real* get_b_x_ptr() { return m_b_x.data(); }
    Real* get_b_x_half_ptr() { return m_b_x_half.data(); }
    Real* get_b_y_ptr() { return m_b_y.data(); }
    Real* get_b_y_half_ptr() { return m_b_y_half.data(); }
    
    // ADE-PML RK4 specific accessors (2D arrays stored as 1D)
    Real* get_a_x_rk4_ptr() { return m_a_x_rk4.data(); }
    Real* get_a_x_half_rk4_ptr() { return m_a_x_half_rk4.data(); }
    Real* get_a_y_rk4_ptr() { return m_a_y_rk4.data(); }
    Real* get_a_y_half_rk4_ptr() { return m_a_y_half_rk4.data(); }
    
    Real* get_b_x_rk4_ptr() { return m_b_x_rk4.data(); }
    Real* get_b_x_half_rk4_ptr() { return m_b_x_half_rk4.data(); }
    Real* get_b_y_rk4_ptr() { return m_b_y_rk4.data(); }
    Real* get_b_y_half_rk4_ptr() { return m_b_y_half_rk4.data(); }
    
    Real* get_memory_dvx_dx_ptr() { return m_memory_dvx_dx.data(); }
    Real* get_memory_dvx_dy_ptr() { return m_memory_dvx_dy.data(); }
    Real* get_memory_dvy_dx_ptr() { return m_memory_dvy_dx.data(); }
    Real* get_memory_dvy_dy_ptr() { return m_memory_dvy_dy.data(); }
    Real* get_memory_dsigmaxx_dx_ptr() { return m_memory_dsigmaxx_dx.data(); }
    Real* get_memory_dsigmayy_dy_ptr() { return m_memory_dsigmayy_dy.data(); }
    Real* get_memory_dsigmaxy_dx_ptr() { return m_memory_dsigmaxy_dx.data(); }
    Real* get_memory_dsigmaxy_dy_ptr() { return m_memory_dsigmaxy_dy.data(); }

    // ADE-PML RK4 specific accessors
    Real* get_lambda_ptr() { return m_lambda.data(); }
    Real* get_mu_ptr() { return m_mu.data(); }
    
    // RK4 intermediate arrays accessors
    Real* get_dvx_ptr() { return m_dvx.data(); }
    Real* get_dvy_ptr() { return m_dvy.data(); }
    Real* get_dsigmaxx_ptr() { return m_dsigmaxx.data(); }
    Real* get_dsigmayy_ptr() { return m_dsigmayy.data(); }
    Real* get_dsigmaxy_ptr() { return m_dsigmaxy.data(); }

    /// Backup current state for RK4 algorithm
    void backup_rk4_state();

    /// Initialize all arrays to zero
    void initialize_fields();

private:
    // Grid dimensions
    Index m_nx, m_ny;
    Real m_deltax, m_deltay;
    
    // Extended grid dimensions for ADE-PML (includes ghost points)
    Index m_nx_extended, m_ny_extended;
    bool m_use_extended_grid;  // Flag to determine which grid to use
    
    // Wave field arrays (stored as 1D vectors in column-major order)
    std::vector<Real> m_vx, m_vy;
    std::vector<Real> m_sigmaxx, m_sigmayy, m_sigmaxy;
    
    // Material property arrays
    std::vector<Real> m_vp, m_vs, m_density;
    std::vector<Real> m_c11, m_c12, m_c22, m_c33;
    std::vector<Real> m_lambda, m_mu;  // ADE-PML requires Lame parameters
    
    // RK4 intermediate arrays (3D: [3 components][nx][ny] stored as 1D)
    std::vector<Real> m_dvx, m_dvy;
    std::vector<Real> m_dsigmaxx, m_dsigmayy, m_dsigmaxy;
    
    // PML profile arrays (1D) - for CPML
    std::vector<Real> m_d_x, m_d_x_half, m_d_y, m_d_y_half;
    std::vector<Real> m_k_x, m_k_x_half, m_k_y, m_k_y_half;
    std::vector<Real> m_alpha_x, m_alpha_x_half, m_alpha_y, m_alpha_y_half;
    std::vector<Real> m_a_x, m_a_x_half, m_a_y, m_a_y_half;
    std::vector<Real> m_b_x, m_b_x_half, m_b_y, m_b_y_half;
    
    // ADE-PML RK4 coefficient arrays (2D: [4 substeps][spatial points])
    std::vector<Real> m_a_x_rk4, m_a_x_half_rk4, m_a_y_rk4, m_a_y_half_rk4;
    std::vector<Real> m_b_x_rk4, m_b_x_half_rk4, m_b_y_rk4, m_b_y_half_rk4;
    
    // PML memory variables (2D)
    std::vector<Real> m_memory_dvx_dx, m_memory_dvx_dy;
    std::vector<Real> m_memory_dvy_dx, m_memory_dvy_dy;
    std::vector<Real> m_memory_dsigmaxx_dx, m_memory_dsigmayy_dy;
    std::vector<Real> m_memory_dsigmaxy_dx, m_memory_dsigmaxy_dy;
};

} // namespace seismic
