#!/usr/bin/env python3
"""
快速查看OpenWave2D结果的简化脚本
"""

import numpy as np
import matplotlib.pyplot as plt
import struct
import glob
import os

def read_snapshot(filename, nx=241, ny=241):
    """读取二进制快照文件"""
    with open(filename, 'rb') as f:
        data = struct.unpack(f'd' * (nx * ny), f.read())
        return np.array(data).reshape((ny, nx))

def quick_view():
    """快速查看最新的结果"""
    output_dir = "build/output"
    
    # 检查输出目录
    if not os.path.exists(output_dir):
        print(f"输出目录不存在: {output_dir}")
        return
    
    # 参数设置
    nx, ny = 241, 241
    deltax, deltay = 10.0, 10.0
    dt = 0.003
    
    print("=" * 50)
    print("OpenWave2D 结果快速查看")
    print("=" * 50)
    
    # 1. 检查地震记录
    seismo_files = glob.glob(f"{output_dir}/seismograms_*.ascii")
    if seismo_files:
        print(f"\n找到 {len(seismo_files)} 个地震记录文件:")
        for f in seismo_files:
            try:
                data = np.loadtxt(f)
                print(f"  {os.path.basename(f)}: {data.shape} (时间步 x 接收器)")
                
                # 绘制第一个接收器的地震记录
                if 'vx' in f:
                    plt.figure(figsize=(12, 4))
                    time_axis = np.arange(data.shape[0]) * dt
                    plt.plot(time_axis, data[:, 0], 'b-', linewidth=1, label='接收器 1')
                    if data.shape[1] > 1:
                        plt.plot(time_axis, data[:, 1], 'r-', linewidth=1, label='接收器 2')
                    plt.xlabel('时间 (s)')
                    plt.ylabel('速度 Vx (m/s)')
                    plt.title('地震记录 - Vx分量')
                    plt.grid(True, alpha=0.3)
                    plt.legend()
                    plt.tight_layout()
                    plt.show()
                    
            except Exception as e:
                print(f"    错误读取 {f}: {e}")
    else:
        print("\n未找到地震记录文件")
    
    # 2. 检查快照文件
    snapshot_files = glob.glob(f"{output_dir}/snapshot_vx_*.bin")
    if snapshot_files:
        # 按时间步排序
        snapshot_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
        print(f"\n找到 {len(snapshot_files)} 个Vx快照文件")
        
        # 显示几个关键时刻的快照
        key_files = [snapshot_files[0], snapshot_files[len(snapshot_files)//4], 
                    snapshot_files[len(snapshot_files)//2], snapshot_files[-1]]
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        x = np.arange(nx) * deltax
        y = np.arange(ny) * deltay
        X, Y = np.meshgrid(x, y)
        
        for i, filename in enumerate(key_files):
            if i >= 4:
                break
                
            try:
                snapshot = read_snapshot(filename, nx, ny)
                timestep = int(filename.split('_')[-1].split('.')[0])
                time = timestep * dt
                
                # 找到颜色范围
                vmax = np.max(np.abs(snapshot))
                if vmax == 0:
                    vmax = 1e-10
                
                im = axes[i].contourf(X, Y, snapshot, levels=20, 
                                    cmap='RdBu_r', vmin=-vmax, vmax=vmax)
                axes[i].set_xlabel('X (m)')
                axes[i].set_ylabel('Y (m)')
                axes[i].set_title(f'Vx - t={time:.3f}s')
                axes[i].set_aspect('equal')
                
                # 标记源位置
                source_x = (nx - 1) * deltax / 2
                source_y = (ny - 1) * deltay / 2
                axes[i].plot(source_x, source_y, 'k*', markersize=8)
                
                plt.colorbar(im, ax=axes[i], shrink=0.8)
                
                # 显示统计信息
                print(f"  时间步 {timestep} (t={time:.3f}s): 最大值={vmax:.2e}")
                
            except Exception as e:
                print(f"    错误读取 {filename}: {e}")
        
        plt.tight_layout()
        plt.show()
        
        # 3. 能量分析
        print("\n计算波场能量...")
        try:
            vx_files = glob.glob(f"{output_dir}/snapshot_vx_*.bin")
            vy_files = glob.glob(f"{output_dir}/snapshot_vy_*.bin")
            
            if len(vx_files) > 0 and len(vy_files) > 0:
                vx_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
                vy_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
                
                times = []
                energies = []
                max_velocities = []
                
                for vx_file, vy_file in zip(vx_files, vy_files):
                    vx_data = read_snapshot(vx_file, nx, ny)
                    vy_data = read_snapshot(vy_file, nx, ny)
                    
                    # 计算总能量 (假设密度=2000)
                    rho = 2000.0
                    energy = 0.5 * rho * np.sum(vx_data**2 + vy_data**2) * deltax * deltay
                    max_vel = np.sqrt(np.max(vx_data**2 + vy_data**2))
                    
                    timestep = int(vx_file.split('_')[-1].split('.')[0])
                    time = timestep * dt
                    
                    times.append(time)
                    energies.append(energy)
                    max_velocities.append(max_vel)
                
                # 绘制能量演化
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
                
                ax1.plot(times, energies, 'b-', linewidth=2)
                ax1.set_xlabel('时间 (s)')
                ax1.set_ylabel('总动能 (J)')
                ax1.set_title('波场总动能演化')
                ax1.grid(True, alpha=0.3)
                
                ax2.plot(times, max_velocities, 'r-', linewidth=2)
                ax2.set_xlabel('时间 (s)')
                ax2.set_ylabel('最大速度 (m/s)')
                ax2.set_title('波场最大速度演化')
                ax2.grid(True, alpha=0.3)
                
                plt.tight_layout()
                plt.show()
                
                print(f"  总动能范围: {min(energies):.2e} - {max(energies):.2e} J")
                print(f"  最大速度范围: {min(max_velocities):.2e} - {max(max_velocities):.2e} m/s")
                
        except Exception as e:
            print(f"  能量分析出错: {e}")
            
    else:
        print("\n未找到快照文件")
    
    print("\n" + "=" * 50)
    print("查看完成!")
    print("=" * 50)

if __name__ == '__main__':
    quick_view()
