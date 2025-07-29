#!/bin/bash

# OpenWave2D 构建脚本
# 用法: ./run_build.sh [clean|rebuild|help|debug]
# 适配Linux环境，支持多种构建工具和编译器

set -e  # 遇到错误时退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 项目根目录
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"

# 系统检测 - Linux专用
detect_system() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        SYSTEM="Linux"
        # 检测具体发行版
        if [ -f /etc/os-release ]; then
            . /etc/os-release
            DISTRO="$NAME"
        elif [ -f /etc/redhat-release ]; then
            DISTRO=$(cat /etc/redhat-release)
        elif [ -f /etc/debian_version ]; then
            DISTRO="Debian $(cat /etc/debian_version)"
        else
            DISTRO="Unknown Linux"
        fi
    else
        echo -e "${RED}错误: 此脚本仅支持Linux系统${NC}"
        echo "当前系统类型: $OSTYPE"
        exit 1
    fi
}

# 检测构建工具
detect_build_tool() {
    if command -v ninja &> /dev/null; then
        BUILD_TOOL="ninja"
        CMAKE_GENERATOR="-G Ninja"
    elif command -v make &> /dev/null; then
        BUILD_TOOL="make"
        CMAKE_GENERATOR=""
    else
        echo -e "${RED}错误: 未找到构建工具 (make 或 ninja)${NC}"
        exit 1
    fi
}

# 检测编译器 - Linux优化版
detect_compiler() {
    # 检测C++编译器（按优先级排序）
    if [[ -n "$CXX" ]]; then
        CXX_COMPILER="$CXX"
        echo -e "${YELLOW}使用环境变量指定的C++编译器: $CXX_COMPILER${NC}"
    elif command -v g++-11 &> /dev/null; then
        CXX_COMPILER="g++-11"
    elif command -v g++-10 &> /dev/null; then
        CXX_COMPILER="g++-10"
    elif command -v g++ &> /dev/null; then
        CXX_COMPILER="g++"
    elif command -v clang++ &> /dev/null; then
        CXX_COMPILER="clang++"
    else
        echo -e "${RED}错误: 未找到C++编译器${NC}"
        echo "请手动安装 g++ 或 clang++"
        exit 1
    fi
    
    # 检测Fortran编译器
    if [[ -n "$FC" ]]; then
        FORTRAN_COMPILER="$FC"
        echo -e "${YELLOW}使用环境变量指定的Fortran编译器: $FORTRAN_COMPILER${NC}"
    elif command -v gfortran-11 &> /dev/null; then
        FORTRAN_COMPILER="gfortran-11"
    elif command -v gfortran-10 &> /dev/null; then
        FORTRAN_COMPILER="gfortran-10"
    elif command -v gfortran &> /dev/null; then
        FORTRAN_COMPILER="gfortran"
    elif command -v ifort &> /dev/null; then
        FORTRAN_COMPILER="ifort"
    else
        echo -e "${RED}错误: 未找到Fortran编译器${NC}"
        echo "请手动安装 gfortran"
        exit 1
    fi
    
    # 显示编译器版本信息
    echo -e "${BLUE}编译器信息:${NC}"
    $CXX_COMPILER --version | head -1
    $FORTRAN_COMPILER --version | head -1
}

# 帮助函数
show_help() {
    echo -e "${BLUE}OpenWave2D 构建脚本 (Linux优化版)${NC}"
    echo ""
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  clean     清理构建目录"
    echo "  rebuild   清理并重新构建"
    echo "  debug     构建Debug版本"
    echo "  release   构建Release版本 (默认)"
    echo "  test      运行测试构建"
    echo "  help      显示此帮助信息"
    echo "  (无参数)   执行标准构建"
    echo ""
    echo "环境变量:"
    echo "  CC        指定C编译器"
    echo "  CXX       指定C++编译器"
    echo "  FC        指定Fortran编译器"
    echo "  BUILD_JOBS 指定并行构建任务数"
    echo ""
    echo "示例:"
    echo "  $0                    # 标准构建"
    echo "  $0 clean              # 清理构建文件"
    echo "  $0 rebuild            # 完全重新构建"
    echo "  $0 debug              # Debug模式构建"
    echo "  CC=clang CXX=clang++ FC=gfortran $0  # 指定编译器"
    echo "  BUILD_JOBS=8 $0       # 使用8个并行任务"
}

# 清理函数
clean_build() {
    echo -e "${YELLOW}清理构建目录...${NC}"
    if [ -d "$BUILD_DIR" ]; then
        rm -rf "$BUILD_DIR"
        echo -e "${GREEN}构建目录已清理${NC}"
    else
        echo -e "${YELLOW}构建目录不存在，无需清理${NC}"
    fi
}

# 创建构建目录
create_build_dir() {
    if [ ! -d "$BUILD_DIR" ]; then
        echo -e "${YELLOW}创建构建目录...${NC}"
        mkdir -p "$BUILD_DIR"
    fi
}

# 配置项目
configure_project() {
    local build_type="${1:-Release}"
    echo -e "${BLUE}配置CMake项目 (${build_type}模式)...${NC}"
    cd "$BUILD_DIR"
    
    # 检查是否有CMake
    if ! command -v cmake &> /dev/null; then
        echo -e "${RED}错误: 未找到CMake，请先安装CMake${NC}"
        exit 1
    fi
    
    # 显示系统信息
    echo -e "${YELLOW}系统信息:${NC}"
    echo "  操作系统: $SYSTEM"
    echo "  发行版: $DISTRO"
    echo "  构建工具: $BUILD_TOOL"
    echo "  C++编译器: $CXX_COMPILER"
    echo "  Fortran编译器: $FORTRAN_COMPILER"
    echo "  构建类型: $build_type"
    echo ""
    
    # 准备CMake参数
    CMAKE_ARGS="-DCMAKE_BUILD_TYPE=$build_type"
    
    # 设置编译器
    if [[ -n "$CXX_COMPILER" ]]; then
        CMAKE_ARGS="$CMAKE_ARGS -DCMAKE_CXX_COMPILER=$CXX_COMPILER"
    fi
    if [[ -n "$FORTRAN_COMPILER" ]]; then
        CMAKE_ARGS="$CMAKE_ARGS -DCMAKE_Fortran_COMPILER=$FORTRAN_COMPILER"
    fi
    
    # 添加生成器
    if [[ -n "$CMAKE_GENERATOR" ]]; then
        CMAKE_ARGS="$CMAKE_ARGS $CMAKE_GENERATOR"
    fi
    
    # 配置项目
    echo -e "${YELLOW}执行: cmake $CMAKE_ARGS $PROJECT_ROOT${NC}"
    cmake $CMAKE_ARGS "$PROJECT_ROOT"
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}CMake配置失败${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}CMake配置完成${NC}"
}

# 构建项目 - Linux优化版
build_project() {
    echo -e "${BLUE}构建项目...${NC}"
    cd "$BUILD_DIR"
    
    # 检测可用的CPU核心数
    if [[ -n "$BUILD_JOBS" ]]; then
        CORES="$BUILD_JOBS"
        echo -e "${YELLOW}使用环境变量指定的 $CORES 个并行任务${NC}"
    elif command -v nproc &> /dev/null; then
        CORES=$(nproc)
        echo -e "${YELLOW}检测到 $CORES 个CPU核心${NC}"
    elif [ -f /proc/cpuinfo ]; then
        CORES=$(grep -c ^processor /proc/cpuinfo)
        echo -e "${YELLOW}从/proc/cpuinfo检测到 $CORES 个CPU核心${NC}"
    else
        CORES=2  # 默认值
        echo -e "${YELLOW}无法检测CPU核心数，使用默认值 $CORES${NC}"
    fi
    
    # 内存检查 - 避免编译时内存不足
    if [ -f /proc/meminfo ]; then
        TOTAL_MEM_KB=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        TOTAL_MEM_GB=$((TOTAL_MEM_KB / 1024 / 1024))
        echo -e "${YELLOW}系统内存: ${TOTAL_MEM_GB}GB${NC}"
        
        # 如果内存小于4GB且CPU核心数大于2，限制并行数
        if [ $TOTAL_MEM_GB -lt 4 ] && [ $CORES -gt 2 ]; then
            CORES=2
            echo -e "${YELLOW}内存不足，限制并行任务数为 $CORES${NC}"
        fi
    fi
    
    # 根据构建工具执行构建
    echo -e "${YELLOW}开始构建 (使用$CORES个并行任务)...${NC}"
    
    if [[ "$BUILD_TOOL" == "ninja" ]]; then
        ninja -j$CORES
    else
        make -j$CORES
    fi
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}构建失败${NC}"
        echo -e "${YELLOW}提示: 如果是内存不足，可以减少并行任务数: BUILD_JOBS=1 $0${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}构建完成！${NC}"
}

# 显示构建结果
show_results() {
    echo -e "${BLUE}构建结果:${NC}"
    echo "  可执行文件: $BUILD_DIR/OpenWave2D"
    echo "  静态库: $BUILD_DIR/libseismic_kernels.a"
    echo "  配置文件: $BUILD_DIR/elastic_homogeneous.json"
    echo "  输出目录: $BUILD_DIR/output/"
    echo "  地震记录目录: $BUILD_DIR/seismograms/"
    echo "  快照目录: $BUILD_DIR/snapshots/"
    echo ""
    echo -e "${GREEN}运行示例:${NC}"
    echo "  cd $BUILD_DIR"
    echo "  ./OpenWave2D elastic_homogeneous.json"
    echo ""
    echo -e "${BLUE}系统要求检查:${NC}"
    
    # 检查依赖库
    if command -v ldd &> /dev/null && [ -f "$BUILD_DIR/OpenWave2D" ]; then
        echo -e "${YELLOW}动态库依赖:${NC}"
        ldd "$BUILD_DIR/OpenWave2D" | grep -E "(not found|libgfortran|libquadmath)" || echo "  所有依赖库已满足"
    fi
}

# 运行快速测试 - Linux增强版
run_test() {
    echo -e "${BLUE}运行快速测试...${NC}"
    cd "$BUILD_DIR"
    
    if [ ! -f "./OpenWave2D" ]; then
        echo -e "${RED}错误: 可执行文件不存在，请先构建项目${NC}"
        exit 1
    fi
    
    # 检查可执行文件权限
    if [ ! -x "./OpenWave2D" ]; then
        echo -e "${YELLOW}添加执行权限...${NC}"
        chmod +x "./OpenWave2D"
    fi
    
    # 检查配置文件
    if [ ! -f "./elastic_homogeneous.json" ]; then
        echo -e "${YELLOW}警告: 配置文件不存在${NC}"
    fi
    
    # 创建测试输出目录
    mkdir -p test_output test_seismograms test_snapshots
    
    # 运行基本测试
    echo -e "${YELLOW}测试可执行文件...${NC}"
    
    # 测试1: 检查可执行文件基本信息
    echo "1. 检查文件信息:"
    file "./OpenWave2D"
    echo ""
    
    # 测试2: 检查动态库依赖
    echo "2. 检查动态库依赖:"
    if command -v ldd &> /dev/null; then
        ldd "./OpenWave2D" | grep -E "(not found|=>)" || echo "  依赖库检查完成"
    fi
    echo ""
    
    # 测试3: 运行帮助信息或简单测试
    echo "3. 运行功能测试:"
    if timeout 10 ./OpenWave2D --help 2>/dev/null; then
        echo -e "${GREEN}帮助信息测试通过${NC}"
    elif timeout 10 ./OpenWave2D -h 2>/dev/null; then
        echo -e "${GREEN}帮助参数测试通过${NC}"
    else
        echo -e "${YELLOW}尝试运行基本测试...${NC}"
        # 如果有配置文件，尝试运行很短时间
        if [ -f "./elastic_homogeneous.json" ]; then
            echo "使用配置文件进行快速测试..."
            timeout 5 ./OpenWave2D elastic_homogeneous.json 2>&1 | head -10
        else
            echo "运行可执行文件基本检查..."
            timeout 3 ./OpenWave2D 2>&1 | head -5
        fi
    fi
    
    echo -e "${GREEN}测试完成${NC}"
    echo ""
    echo -e "${BLUE}如需完整测试，请运行:${NC}"
    echo "  cd $BUILD_DIR"
    echo "  ./OpenWave2D elastic_homogeneous.json"
}

# 主逻辑
main() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}    OpenWave2D Linux构建系统${NC}"
    echo -e "${BLUE}========================================${NC}"
    
    # 系统检测
    detect_system
    detect_build_tool
    detect_compiler
    
    case "${1:-}" in
        "clean")
            clean_build
            ;;
        "rebuild")
            clean_build
            create_build_dir
            configure_project "Release"
            build_project
            show_results
            ;;
        "debug")
            create_build_dir
            configure_project "Debug"
            build_project
            show_results
            ;;
        "release")
            create_build_dir
            configure_project "Release"
            build_project
            show_results
            ;;
        "test")
            if [ ! -f "$BUILD_DIR/OpenWave2D" ]; then
                create_build_dir
                configure_project "Release"
                build_project
            fi
            run_test
            ;;
        "help"|"-h"|"--help")
            show_help
            ;;
        "")
            create_build_dir
            configure_project "Release"
            build_project
            show_results
            ;;
        *)
            echo -e "${RED}未知选项: $1${NC}"
            echo ""
            show_help
            exit 1
            ;;
    esac
}

# 检查是否在项目根目录
if [ ! -f "$PROJECT_ROOT/CMakeLists.txt" ]; then
    echo -e "${RED}错误: 未在项目根目录中找到CMakeLists.txt${NC}"
    echo "请确保在OpenWave2D项目根目录中运行此脚本"
    exit 1
fi

# 执行主函数
main "$@"
