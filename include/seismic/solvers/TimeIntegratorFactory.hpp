#pragma once

#include <memory>
#include <string>
#include "seismic/core/Config.hpp"
#include "TimeIntegrator.hpp"

namespace seismic {

/**
 * @class TimeIntegratorFactory
 * @brief 时间积分器工厂类
 * 
 * 这个工厂类负责根据配置参数创建适当的时间积分器实例。
 * 通过工厂模式，我们可以轻松地添加新的积分器类型，而无需
 * 修改主调度器的代码。
 */
class TimeIntegratorFactory {
public:
    /**
     * @brief 根据PML类型和时间积分方案创建时间积分器实例
     * @param time_config 时间配置参数
     * @param pml_config PML配置参数
     * @return 指向TimeIntegrator实例的唯一指针
     * 
     * 支持的组合：
     * - PML类型 "cpml" + 任何时间方案 -> CPMLIntegrator
     * - PML类型 "adepml" + 时间方案 "rk4" -> RK4Integrator
     */
    static std::unique_ptr<TimeIntegrator> create(const TimeConfig& time_config,
                                                 const PMLConfig& pml_config);

    /**
     * @brief 检查给定的配置组合是否支持
     * @param time_config 时间配置参数
     * @param pml_config PML配置参数
     * @return 如果支持返回true，否则返回false
     */
    static bool is_supported(const TimeConfig& time_config,
                            const PMLConfig& pml_config);

    /**
     * @brief 获取所有支持的积分器类型列表（用于错误信息和帮助）
     * @return 支持的积分器类型字符串
     */
    static std::string get_supported_types();

private:
    // 工厂类不需要实例化
    TimeIntegratorFactory() = delete;
};

/**
 * @brief 便利的工厂函数（为了保持向后兼容性）
 * @param time_config 时间配置
 * @param pml_config PML配置
 * @return 指向TimeIntegrator实例的唯一指针
 */
std::unique_ptr<TimeIntegrator> create_time_integrator(const TimeConfig& time_config,
                                                      const PMLConfig& pml_config);

} // namespace seismic
