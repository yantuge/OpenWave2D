#pragma once

#include "seismic/core/Grid.hpp"
#include "seismic/core/Config.hpp"
#include "seismic/sources/SourceTimeFunction.hpp"

// Fortran内核函数声明，避免在头文件中包含Dispatcher.hpp
extern "C" {
    void cpml_anisotropic_step_c(int nx, int ny, double deltax, double deltay, double dt, int it,
                                int isource, int jsource, double force_x, double force_y,
                                double* c11, double* c12, double* c22, double* c33, double* rho,
                                double* vx, double* vy, double* sigmaxx, double* sigmayy, double* sigmaxy,
                                double* d_x, double* d_x_half, double* d_y, double* d_y_half,
                                double* k_x, double* k_x_half, double* k_y, double* k_y_half,
                                double* alpha_x, double* alpha_x_half, double* alpha_y, double* alpha_y_half,
                                double* a_x, double* a_x_half, double* a_y, double* a_y_half,
                                double* b_x, double* b_x_half, double* b_y, double* b_y_half,
                                double* memory_dvx_dx, double* memory_dvx_dy,
                                double* memory_dvy_dx, double* memory_dvy_dy,
                                double* memory_dsigmaxx_dx, double* memory_dsigmayy_dy,
                                double* memory_dsigmaxy_dx, double* memory_dsigmaxy_dy);
    
    void adepml_rk4_step_c(int nx, int ny, double deltax, double deltay, double dt, int rk_substep,
                          int isource, int jsource, double force_x, double force_y,
                          double* lambda, double* mu, double* rho,
                          double* vx, double* vy, double* sigmaxx, double* sigmayy, double* sigmaxy,
                          double* dvx, double* dvy, double* dsigmaxx, double* dsigmayy, double* dsigmaxy,
                          double* d_x, double* d_x_half, double* d_y, double* d_y_half,
                          double* k_x, double* k_x_half, double* k_y, double* k_y_half,
                          double* alpha_x, double* alpha_x_half, double* alpha_y, double* alpha_y_half,
                          double* a_x_rk4, double* a_x_half_rk4, double* a_y_rk4, double* a_y_half_rk4,
                          double* b_x_rk4, double* b_x_half_rk4, double* b_y_rk4, double* b_y_half_rk4,
                          double* memory_dvx_dx, double* memory_dvx_dy,
                          double* memory_dvy_dx, double* memory_dvy_dy,
                          double* memory_dsigmaxx_dx, double* memory_dsigmayy_dy,
                          double* memory_dsigmaxy_dx, double* memory_dsigmaxy_dy);
}

namespace seismic {

/**
 * @class TimeIntegrator
 * @brief 时间积分方案的抽象基类 (策略接口)
 * 
 * 这个接口定义了所有时间积分方案必须遵守的契约。通过策略模式，
 * 我们可以在运行时根据配置选择不同的时间积分算法，而无需修改
 * 主调度器的代码。
 */
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    /**
     * @brief 执行单个时间步的积分
     * @param grid 网格数据对象，包含所有场变量和PML参数
     * @param source_config 震源配置参数
     * @param pml_config PML边界条件配置
     * @param dt 时间步长
     * @param it 当前时间步索引
     * 
     * 每个积分器实现必须：
     * 1. 计算当前时间步的震源力
     * 2. 调用适当的Fortran内核进行波场更新
     * 3. 处理任何特定于该积分方案的逻辑
     */
    virtual void step(Grid& grid, 
                     const SourceConfig& source_config,
                     const PMLConfig& pml_config, 
                     Real dt, 
                     Index it) = 0;

    /**
     * @brief 获取积分器的名称（用于日志和调试）
     */
    virtual std::string get_name() const = 0;
    
    /**
     * @brief 获取积分器的稳定性条件系数
     * 不同的时间积分方案可能有不同的稳定性要求
     */
    virtual Real get_stability_factor() const = 0;
};

} // namespace seismic
