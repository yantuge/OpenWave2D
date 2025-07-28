#pragma once

#include "TimeIntegrator.hpp"

namespace seismic {

/**
 * @class CPMLIntegrator
 * @brief CPML (Convolutional Perfectly Matched Layer) 显式时间积分器
 * 
 * 这个类实现了CPML的显式时间积分方案。与RK4不同，这是一个
 * 一步显式方法，每个时间步只需要调用一次Fortran内核。
 * 
 * CPML的特点：
 * - 计算效率高（每时间步只需一次内核调用）
 * - 内存需求较低
 * - 适用于各向异性介质
 */
class CPMLIntegrator : public TimeIntegrator {
public:
    /**
     * @brief 构造函数
     */
    CPMLIntegrator() = default;
    
    /**
     * @brief 析构函数
     */
    ~CPMLIntegrator() override = default;

    /**
     * @brief 执行CPML显式时间步积分
     */
    void step(Grid& grid, 
             const SourceConfig& source_config,
             const PMLConfig& pml_config, 
             Real dt, 
             Index it) override;

    /**
     * @brief 获取积分器名称
     */
    std::string get_name() const override { return "CPML"; }
    
    /**
     * @brief 获取CPML的稳定性系数
     * CPML显式方法需要严格满足CFL条件
     */
    Real get_stability_factor() const override { return 1.0; }

private:
    /**
     * @brief 计算震源在指定时间的贡献
     */
    void calculate_source_forces_at_time(Real time,
                                        const SourceConfig& source_config,
                                        Real& force_x, 
                                        Real& force_y) const;
};

} // namespace seismic
