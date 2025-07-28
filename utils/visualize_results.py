#!/usr/bin/env python3
"""
OpenWave2D结果可视化工具

该脚本用于可视化OpenWave2D地震模拟的结果，包括：
1. 地震记录（seismograms）的显示
2. 波场快照（snapshots）的显示
3. 波场传播动画的生成

用法:
    python visualize_results.py --help
    
    查看地震记录:
    python visualize_results.py --seismograms
    
    查看波场快照:
    python visualize_results.py --snapshots
    
    生成动画:
    python visualize_results.py --animation
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse
import os
import glob
import struct
from pathlib import Path

class OpenWave2DVisualizer:
    def __init__(self, output_dir="build/output", config_file="build/elastic_homogeneous.json"):
        self.output_dir = Path(output_dir)
        self.config_file = Path(config_file)
        
        # 默认参数（如果无法读取配置文件）
        self.nx = 241
        self.ny = 241
        self.deltax = 10.0
        self.deltay = 10.0
        self.dt = 0.003
        self.nstep = 2501
        
        # 尝试读取配置文件参数
        self.load_config()
        
    def load_config(self):
        """从JSON配置文件读取模拟参数"""
        try:
            import json
            with open(self.config_file, 'r') as f:
                config = json.load(f)
                
            self.nx = config['grid']['nx']
            self.ny = config['grid']['ny']
            self.deltax = config['grid']['deltax']
            self.deltay = config['grid']['deltay']
            self.dt = config['time']['dt']
            self.nstep = config['time']['nstep']
            
            print(f"从配置文件读取参数: nx={self.nx}, ny={self.ny}, dt={self.dt}")
            
        except Exception as e:
            print(f"警告: 无法读取配置文件 {self.config_file}: {e}")
            print(f"使用默认参数: nx={self.nx}, ny={self.ny}, dt={self.dt}")
    
    def read_seismograms(self):
        """读取地震记录数据"""
        seismograms = {}
        
        # 查找地震记录文件
        for component in ['vx', 'vy']:
            filename = self.output_dir / f"seismograms_{component}.ascii"
            if filename.exists():
                try:
                    # 读取ASCII格式的地震记录
                    data = np.loadtxt(filename)
                    seismograms[component] = data
                    print(f"读取地震记录: {filename}, 形状: {data.shape}")
                except Exception as e:
                    print(f"错误: 无法读取 {filename}: {e}")
            else:
                print(f"警告: 地震记录文件不存在: {filename}")
                
        return seismograms
    
    def read_snapshot(self, filename):
        """读取二进制快照文件"""
        try:
            with open(filename, 'rb') as f:
                # 二进制文件包含 nx*ny 个double值
                data = struct.unpack(f'd' * (self.nx * self.ny), f.read())
                # 重塑为2D数组 (注意：Fortran列主序 vs C行主序)
                snapshot = np.array(data).reshape((self.ny, self.nx))
                return snapshot
        except Exception as e:
            print(f"错误: 无法读取快照文件 {filename}: {e}")
            return None
    
    def get_snapshot_files(self, component='vx'):
        """获取所有快照文件的列表"""
        pattern = str(self.output_dir / f"snapshot_{component}_*.bin")
        files = sorted(glob.glob(pattern), key=lambda x: int(x.split('_')[-1].split('.')[0]))
        return files
    
    def plot_seismograms(self):
        """绘制地震记录"""
        seismograms = self.read_seismograms()
        
        if not seismograms:
            print("没有找到地震记录数据")
            return
            
        fig, axes = plt.subplots(len(seismograms), 1, figsize=(12, 4*len(seismograms)))
        if len(seismograms) == 1:
            axes = [axes]
            
        time_axis = np.arange(seismograms[list(seismograms.keys())[0]].shape[0]) * self.dt
        
        for i, (component, data) in enumerate(seismograms.items()):
            ax = axes[i]
            
            # 绘制所有接收器的地震记录
            for receiver in range(data.shape[1]):
                ax.plot(time_axis, data[:, receiver], 
                       label=f'接收器 {receiver+1}', linewidth=0.8)
            
            ax.set_xlabel('时间 (s)')
            ax.set_ylabel(f'速度 {component.upper()} (m/s)')
            ax.set_title(f'地震记录 - {component.upper()} 分量')
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
        plt.tight_layout()
        plt.show()
    
    def plot_snapshots(self, component='vx', timesteps=None):
        """绘制波场快照"""
        snapshot_files = self.get_snapshot_files(component)
        
        if not snapshot_files:
            print(f"没有找到 {component} 分量的快照文件")
            return
            
        if timesteps is None:
            # 显示前6个快照
            timesteps = snapshot_files[:6]
        else:
            timesteps = [f for f in snapshot_files if any(str(t) in f for t in timesteps)]
            
        n_plots = len(timesteps)
        cols = 3
        rows = (n_plots + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(15, 5*rows))
        if rows == 1:
            axes = axes.reshape(1, -1)
        if cols == 1:
            axes = axes.reshape(-1, 1)
            
        # 坐标轴
        x = np.arange(self.nx) * self.deltax
        y = np.arange(self.ny) * self.deltay
        X, Y = np.meshgrid(x, y)
        
        for i, filename in enumerate(timesteps):
            row, col = divmod(i, cols)
            ax = axes[row, col]
            
            snapshot = self.read_snapshot(filename)
            if snapshot is not None:
                # 提取时间步
                timestep = int(filename.split('_')[-1].split('.')[0])
                time = timestep * self.dt
                
                # 找到合适的颜色范围
                vmax = np.max(np.abs(snapshot))
                if vmax == 0:
                    vmax = 1e-10
                    
                im = ax.contourf(X, Y, snapshot, levels=50, 
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax)
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Y (m)')
                ax.set_title(f'{component.upper()} - t={time:.3f}s (步骤 {timestep})')
                ax.set_aspect('equal')
                
                # 添加颜色条
                plt.colorbar(im, ax=ax, shrink=0.8)
                
                # 标记源位置 (假设在中心)
                source_x = (self.nx - 1) * self.deltax / 2
                source_y = (self.ny - 1) * self.deltay / 2
                ax.plot(source_x, source_y, 'k*', markersize=10, label='源')
        
        # 隐藏多余的子图
        for i in range(n_plots, rows * cols):
            row, col = divmod(i, cols)
            axes[row, col].axis('off')
            
        plt.tight_layout()
        plt.show()
    
    def create_animation(self, component='vx', output_filename='wave_animation.mp4', 
                        interval=100, save_gif=True):
        """创建波场传播动画"""
        snapshot_files = self.get_snapshot_files(component)
        
        if not snapshot_files:
            print(f"没有找到 {component} 分量的快照文件")
            return
            
        print(f"创建动画，共 {len(snapshot_files)} 帧...")
        
        # 读取第一个快照来设置图形
        first_snapshot = self.read_snapshot(snapshot_files[0])
        if first_snapshot is None:
            return
            
        # 找到所有快照的全局最大值来设置颜色范围
        global_max = 0
        for filename in snapshot_files[::5]:  # 采样检查最大值
            snapshot = self.read_snapshot(filename)
            if snapshot is not None:
                global_max = max(global_max, np.max(np.abs(snapshot)))
                
        if global_max == 0:
            global_max = 1e-10
            
        # 设置图形
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # 坐标轴
        x = np.arange(self.nx) * self.deltax
        y = np.arange(self.ny) * self.deltay
        X, Y = np.meshgrid(x, y)
        
        # 初始化绘图
        im = ax.contourf(X, Y, first_snapshot, levels=50, 
                        cmap='RdBu_r', vmin=-global_max, vmax=global_max)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_aspect('equal')
        
        # 添加颜色条
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(f'速度 {component.upper()} (m/s)')
        
        # 标记源位置
        source_x = (self.nx - 1) * self.deltax / 2
        source_y = (self.ny - 1) * self.deltay / 2
        ax.plot(source_x, source_y, 'k*', markersize=15, label='震源')
        ax.legend()
        
        title = ax.set_title('')
        
        def animate(frame):
            if frame < len(snapshot_files):
                filename = snapshot_files[frame]
                snapshot = self.read_snapshot(filename)
                
                if snapshot is not None:
                    # 清除之前的等高线
                    ax.clear()
                    
                    # 重新绘制
                    im = ax.contourf(X, Y, snapshot, levels=50, 
                                   cmap='RdBu_r', vmin=-global_max, vmax=global_max)
                    
                    # 提取时间步
                    timestep = int(filename.split('_')[-1].split('.')[0])
                    time = timestep * self.dt
                    
                    ax.set_xlabel('X (m)')
                    ax.set_ylabel('Y (m)')
                    ax.set_title(f'波场传播 ({component.upper()}) - t={time:.3f}s')
                    ax.set_aspect('equal')
                    
                    # 重新添加源标记
                    ax.plot(source_x, source_y, 'k*', markersize=15)
                    
            return []
        
        # 创建动画
        anim = animation.FuncAnimation(fig, animate, frames=len(snapshot_files),
                                     interval=interval, blit=False, repeat=True)
        
        # 保存动画
        if save_gif:
            gif_filename = output_filename.replace('.mp4', '.gif')
            print(f"保存GIF动画: {gif_filename}")
            anim.save(gif_filename, writer='pillow', fps=10)
            
        try:
            print(f"保存MP4动画: {output_filename}")
            anim.save(output_filename, writer='ffmpeg', fps=10)
        except Exception as e:
            print(f"无法保存MP4文件 (需要ffmpeg): {e}")
            
        plt.show()
        return anim
    
    def analyze_energy(self):
        """分析波场能量演化"""
        snapshot_files_vx = self.get_snapshot_files('vx')
        snapshot_files_vy = self.get_snapshot_files('vy')
        
        if not snapshot_files_vx or not snapshot_files_vy:
            print("需要vx和vy快照文件来计算能量")
            return
            
        times = []
        kinetic_energy = []
        
        for i, (file_vx, file_vy) in enumerate(zip(snapshot_files_vx, snapshot_files_vy)):
            vx_data = self.read_snapshot(file_vx)
            vy_data = self.read_snapshot(file_vy)
            
            if vx_data is not None and vy_data is not None:
                # 计算动能密度 (假设密度为2000 kg/m³)
                rho = 2000.0
                ke = 0.5 * rho * (vx_data**2 + vy_data**2)
                total_ke = np.sum(ke) * self.deltax * self.deltay
                
                timestep = int(file_vx.split('_')[-1].split('.')[0])
                time = timestep * self.dt
                
                times.append(time)
                kinetic_energy.append(total_ke)
        
        # 绘制能量演化
        plt.figure(figsize=(10, 6))
        plt.plot(times, kinetic_energy, 'b-', linewidth=2)
        plt.xlabel('时间 (s)')
        plt.ylabel('总动能 (J)')
        plt.title('波场总动能随时间演化')
        plt.grid(True, alpha=0.3)
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='OpenWave2D结果可视化工具')
    parser.add_argument('--output-dir', default='build/output', 
                       help='输出目录路径 (默认: build/output)')
    parser.add_argument('--config', default='build/elastic_homogeneous.json',
                       help='配置文件路径 (默认: build/elastic_homogeneous.json)')
    parser.add_argument('--seismograms', action='store_true',
                       help='显示地震记录')
    parser.add_argument('--snapshots', action='store_true',
                       help='显示波场快照')
    parser.add_argument('--animation', action='store_true',
                       help='创建波场传播动画')
    parser.add_argument('--energy', action='store_true',
                       help='分析波场能量演化')
    parser.add_argument('--component', default='vx', choices=['vx', 'vy'],
                       help='选择要可视化的速度分量 (默认: vx)')
    parser.add_argument('--output-animation', default='wave_animation.mp4',
                       help='动画输出文件名 (默认: wave_animation.mp4)')
    
    args = parser.parse_args()
    
    # 创建可视化器
    visualizer = OpenWave2DVisualizer(args.output_dir, args.config)
    
    # 如果没有指定任何选项，显示帮助
    if not any([args.seismograms, args.snapshots, args.animation, args.energy]):
        print("请选择要执行的操作:")
        print("  --seismograms  : 显示地震记录")
        print("  --snapshots    : 显示波场快照")
        print("  --animation    : 创建波场传播动画")
        print("  --energy       : 分析波场能量演化")
        print("\n示例:")
        print("  python visualize_results.py --seismograms")
        print("  python visualize_results.py --snapshots --component vx")
        print("  python visualize_results.py --animation --component vy")
        return
    
    # 执行用户选择的操作
    if args.seismograms:
        print("显示地震记录...")
        visualizer.plot_seismograms()
        
    if args.snapshots:
        print(f"显示波场快照 ({args.component} 分量)...")
        visualizer.plot_snapshots(component=args.component)
        
    if args.animation:
        print(f"创建波场传播动画 ({args.component} 分量)...")
        visualizer.create_animation(component=args.component, 
                                   output_filename=args.output_animation)
        
    if args.energy:
        print("分析波场能量演化...")
        visualizer.analyze_energy()

if __name__ == '__main__':
    main()
