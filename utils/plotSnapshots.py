import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import glob
from pathlib import Path

def read_su_snapshot(file_path):
    """
    读取SU格式的快照文件
    
    Parameters:
    file_path: SU文件路径
    
    Returns:
    data: 2D numpy数组，形状为(ny, nx)
    """
    try:
        # 获取文件信息
        file_size = os.path.getsize(file_path)
        
        # SU格式参数（根据我们的C++实现）
        header_size = 240  # SU头大小（字节）
        nx, ny = 241, 241  # 网格尺寸
        bytes_per_sample = 4  # float32
        
        # 计算每道的字节数
        samples_per_trace = ny
        bytes_per_trace = header_size + samples_per_trace * bytes_per_sample
        expected_size = nx * bytes_per_trace
        
        if file_size != expected_size:
            print(f"警告：文件大小({file_size})与期望大小({expected_size})不匹配")
        
        # 读取数据
        data = np.zeros((ny, nx), dtype=np.float32)
        
        with open(file_path, 'rb') as f:
            for i in range(nx):
                # 跳过SU头
                f.seek(i * bytes_per_trace + header_size)
                # 读取一列数据
                column_data = np.frombuffer(f.read(samples_per_trace * bytes_per_sample), dtype=np.float32)
                if len(column_data) == ny:
                    data[:, i] = column_data
                else:
                    print(f"警告：第{i}列数据长度({len(column_data)})不正确")
                    break
        
        return data
        
    except Exception as e:
        print(f"读取文件{file_path}时出错: {e}")
        return None

def plot_single_snapshot(file_path, component='Vx', save_fig=False, output_dir='./plots'):
    """
    绘制单个快照
    
    Parameters:
    file_path: SU文件路径
    component: 分量名称('Vx'或'Vy')
    save_fig: 是否保存图片
    output_dir: 输出目录
    """
    data = read_su_snapshot(file_path)
    if data is None:
        return
    
    # 提取时间步信息
    filename = os.path.basename(file_path)
    time_step = int(filename[:5])
    time = time_step * 0.003  # DELTAT = 0.003秒
    
    # 设置网格参数
    dx, dy = 10.0, 10.0  # 网格间距（米）
    nx, ny = data.shape[1], data.shape[0]
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    
    # 创建图形
    plt.figure(figsize=(12, 10))
    
    # 绘制波场
    vmax = np.max(np.abs(data))
    plt.imshow(data, extent=[0, (nx-1)*dx, (ny-1)*dy, 0], 
               cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
    
    plt.colorbar(label=f'{component} Velocity (m/s)', shrink=0.8)
    plt.title(f'{component} Wave Field - Time Step {time_step} (t = {time:.3f}s)', fontsize=14)
    plt.xlabel('X Distance (m)', fontsize=12)
    plt.ylabel('Y Distance (m)', fontsize=12)
    
    # 添加网格
    plt.grid(True, alpha=0.3)
    
    # 添加震源位置标记（假设在中心）
    source_x, source_y = 1190, 1190  # 根据C++代码中的震源位置
    plt.plot(source_x, source_y, 'k*', markersize=15, label='Source')
    plt.legend()
    
    plt.tight_layout()
    
    if save_fig:
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f'{component}_snapshot_{time_step:05d}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"图片已保存: {output_file}")
    
    plt.show()

def plot_snapshot_comparison(vx_file, vy_file, save_fig=False, output_dir='./plots'):
    """
    并排比较Vx和Vy快照
    
    Parameters:
    vx_file: Vx分量文件路径
    vy_file: Vy分量文件路径
    save_fig: 是否保存图片
    output_dir: 输出目录
    """
    vx_data = read_su_snapshot(vx_file)
    vy_data = read_su_snapshot(vy_file)
    
    if vx_data is None or vy_data is None:
        return
    
    # 提取时间步信息
    filename = os.path.basename(vx_file)
    time_step = int(filename[:5])
    time = time_step * 0.003
    
    # 设置网格参数
    dx, dy = 10.0, 10.0
    nx, ny = vx_data.shape[1], vx_data.shape[0]
    
    # 统一颜色范围
    vmax = max(np.max(np.abs(vx_data)), np.max(np.abs(vy_data)))
    
    # 创建子图
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Vx分量
    im1 = ax1.imshow(vx_data, extent=[0, (nx-1)*dx, (ny-1)*dy, 0], 
                     cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
    ax1.set_title(f'Vx Component (t = {time:.3f}s)')
    ax1.set_xlabel('X Distance (m)')
    ax1.set_ylabel('Y Distance (m)')
    ax1.plot(1190, 1190, 'k*', markersize=10)
    
    # Vy分量
    im2 = ax2.imshow(vy_data, extent=[0, (nx-1)*dx, (ny-1)*dy, 0], 
                     cmap='seismic', vmin=-vmax, vmax=vmax, aspect='equal')
    ax2.set_title(f'Vy Component (t = {time:.3f}s)')
    ax2.set_xlabel('X Distance (m)')
    ax2.set_ylabel('Y Distance (m)')
    ax2.plot(1190, 1190, 'k*', markersize=10)
    
    # 速度幅值
    velocity_magnitude = np.sqrt(vx_data**2 + vy_data**2)
    im3 = ax3.imshow(velocity_magnitude, extent=[0, (nx-1)*dx, (ny-1)*dy, 0], 
                     cmap='hot', vmin=0, vmax=np.max(velocity_magnitude), aspect='equal')
    ax3.set_title(f'Velocity Magnitude (t = {time:.3f}s)')
    ax3.set_xlabel('X Distance (m)')
    ax3.set_ylabel('Y Distance (m)')
    ax3.plot(1190, 1190, 'k*', markersize=10)
    
    # 添加颜色条
    plt.colorbar(im1, ax=ax1, label='Vx (m/s)', shrink=0.8)
    plt.colorbar(im2, ax=ax2, label='Vy (m/s)', shrink=0.8)
    plt.colorbar(im3, ax=ax3, label='|V| (m/s)', shrink=0.8)
    
    plt.suptitle(f'Wave Field Snapshot - Time Step {time_step}', fontsize=16)
    plt.tight_layout()
    
    if save_fig:
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f'comparison_snapshot_{time_step:05d}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"比较图已保存: {output_file}")
    
    plt.show()

def create_animation(snapshots_dir, component='Vx', output_file='wave_animation.gif', 
                    fps=5, save_frames=False, output_dir='./plots'):
    """
    创建波场传播动画
    
    Parameters:
    snapshots_dir: 快照文件目录
    component: 分量('Vx'或'Vy')
    output_file: 输出动画文件名
    fps: 帧率
    save_frames: 是否保存单独的帧
    output_dir: 输出目录
    """
    # 查找所有快照文件
    pattern = os.path.join(snapshots_dir, f'*snap{component}.su')
    snapshot_files = sorted(glob.glob(pattern))
    
    if not snapshot_files:
        print(f"在{snapshots_dir}中没有找到{component}快照文件")
        return
    
    print(f"找到{len(snapshot_files)}个{component}快照文件")
    
    # 读取第一个文件确定数据形状
    first_data = read_su_snapshot(snapshot_files[0])
    if first_data is None:
        return
    
    ny, nx = first_data.shape
    dx, dy = 10.0, 10.0
    
    # 设置图形
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # 计算全局颜色范围
    print("计算全局颜色范围...")
    all_max = 0
    for i, file_path in enumerate(snapshot_files[::5]):  # 采样计算以提高速度
        data = read_su_snapshot(file_path)
        if data is not None:
            all_max = max(all_max, np.max(np.abs(data)))
        if i % 10 == 0:
            print(f"处理进度: {i*5}/{len(snapshot_files)}")
    
    # 初始化图像
    im = ax.imshow(first_data, extent=[0, (nx-1)*dx, (ny-1)*dy, 0], 
                   cmap='seismic', vmin=-all_max, vmax=all_max, aspect='equal')
    
    ax.set_xlabel('X Distance (m)')
    ax.set_ylabel('Y Distance (m)')
    ax.plot(1190, 1190, 'k*', markersize=15, label='Source')
    ax.legend()
    
    cbar = plt.colorbar(im, label=f'{component} Velocity (m/s)')
    title = ax.set_title('')
    
    def animate(frame):
        if frame < len(snapshot_files):
            data = read_su_snapshot(snapshot_files[frame])
            if data is not None:
                im.set_array(data)
                
                # 更新标题
                filename = os.path.basename(snapshot_files[frame])
                time_step = int(filename[:5])
                time = time_step * 0.003
                title.set_text(f'{component} Wave Field - Time Step {time_step} (t = {time:.3f}s)')
                
                if save_frames:
                    frames_dir = os.path.join(output_dir, 'animation_frames')
                    os.makedirs(frames_dir, exist_ok=True)
                    frame_file = os.path.join(frames_dir, f'frame_{frame:04d}.png')
                    plt.savefig(frame_file, dpi=100, bbox_inches='tight')
        
        return [im, title]
    
    # 创建动画
    print("创建动画...")
    anim = animation.FuncAnimation(fig, animate, frames=len(snapshot_files), 
                                 interval=1000//fps, blit=False, repeat=True)
    
    # 保存动画
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)
    print(f"保存动画到 {output_path}...")
    anim.save(output_path, writer='pillow', fps=fps, dpi=100)
    print(f"动画已保存: {output_path}")
    
    plt.show()

def main():
    """
    主函数 - 提供交互式菜单
    """
    snapshots_dir = r'E:\OpenWave2D\build\output\snapshots'
    
    print("=" * 50)
    print("OpenWave++ 快照可视化工具")
    print("=" * 50)
    print("1. 绘制单个Vx快照")
    print("2. 绘制单个Vy快照") 
    print("3. 绘制Vx和Vy对比快照")
    print("4. 创建Vx波场传播动画")
    print("5. 创建Vy波场传播动画")
    print("6. 列出所有可用的快照文件")
    print("0. 退出")
    print("=" * 50)
    
    while True:
        choice = input("\n请选择功能 (0-6): ").strip()
        
        if choice == '0':
            print("退出程序")
            break
            
        elif choice == '1':
            # 列出Vx文件
            vx_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVx.su')))
            if not vx_files:
                print("没有找到Vx快照文件")
                continue
            
            print(f"找到{len(vx_files)}个Vx快照文件")
            for i, f in enumerate(vx_files[:10]):  # 只显示前10个
                time_step = int(os.path.basename(f)[:5])
                print(f"{i}: {os.path.basename(f)} (t = {time_step*0.003:.3f}s)")
            
            if len(vx_files) > 10:
                print(f"... 还有{len(vx_files)-10}个文件")
            
            try:
                idx = int(input(f"选择文件索引 (0-{min(9, len(vx_files)-1)}): "))
                if 0 <= idx < len(vx_files):
                    plot_single_snapshot(vx_files[idx], 'Vx', save_fig=True)
            except ValueError:
                print("无效输入")
                
        elif choice == '2':
            # 类似处理Vy文件
            vy_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVy.su')))
            if not vy_files:
                print("没有找到Vy快照文件")
                continue
                
            print(f"找到{len(vy_files)}个Vy快照文件")
            for i, f in enumerate(vy_files[:10]):
                time_step = int(os.path.basename(f)[:5])
                print(f"{i}: {os.path.basename(f)} (t = {time_step*0.003:.3f}s)")
            
            if len(vy_files) > 10:
                print(f"... 还有{len(vy_files)-10}个文件")
            
            try:
                idx = int(input(f"选择文件索引 (0-{min(9, len(vy_files)-1)}): "))
                if 0 <= idx < len(vy_files):
                    plot_single_snapshot(vy_files[idx], 'Vy', save_fig=True)
            except ValueError:
                print("无效输入")
                
        elif choice == '3':
            # 对比快照
            vx_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVx.su')))
            vy_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVy.su')))
            
            if not vx_files or not vy_files:
                print("没有找到足够的快照文件")
                continue
            
            print(f"找到{len(vx_files)}对快照文件")
            for i, f in enumerate(vx_files[:10]):
                time_step = int(os.path.basename(f)[:5])
                print(f"{i}: 时间步{time_step} (t = {time_step*0.003:.3f}s)")
            
            try:
                idx = int(input(f"选择时间步索引 (0-{min(9, len(vx_files)-1)}): "))
                if 0 <= idx < len(vx_files):
                    plot_snapshot_comparison(vx_files[idx], vy_files[idx], save_fig=True)
            except ValueError:
                print("无效输入")
                
        elif choice == '4':
            # Vx动画
            output_file = input("输入动画文件名 (默认: vx_wave_animation.gif): ").strip()
            if not output_file:
                output_file = 'vx_wave_animation.gif'
            create_animation(snapshots_dir, 'Vx', output_file, output_dir='./plots')
            
        elif choice == '5':
            # Vy动画
            output_file = input("输入动画文件名 (默认: vy_wave_animation.gif): ").strip()
            if not output_file:
                output_file = 'vy_wave_animation.gif'
            create_animation(snapshots_dir, 'Vy', output_file, output_dir='./plots')
            
        elif choice == '6':
            # 列出文件
            vx_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVx.su')))
            vy_files = sorted(glob.glob(os.path.join(snapshots_dir, '*snapVy.su')))
            
            print(f"\n在 {snapshots_dir} 中找到:")
            print(f"- {len(vx_files)} 个 Vx 快照文件")
            print(f"- {len(vy_files)} 个 Vy 快照文件")
            
            if vx_files:
                first_time = int(os.path.basename(vx_files[0])[:5]) * 0.003
                last_time = int(os.path.basename(vx_files[-1])[:5]) * 0.003
                print(f"时间范围: {first_time:.3f}s - {last_time:.3f}s")
        
        else:
            print("无效选择，请重试")

if __name__ == "__main__":
    main()
