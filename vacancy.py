mport sys
import os

def process_frames(input_file, remove_counts):
    # 预处理：读取所有帧到内存（适合中小型文件）
    frames = []
    with open(input_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            n_atoms = int(line.strip())
            comment = f.readline().strip()
            atoms = []
            for _ in range(n_atoms):
                parts = f.readline().split()
                elem, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                atoms.append((elem, x, y, z))
            frames.append((n_atoms, comment, atoms))
    
    # 为每个删除数量创建输出文件
    for count in remove_counts:
        output_file = f"output_{count}.xyz"
        with open(output_file, 'w') as outfile:
            for frame in frames:
                n_atoms_orig, comment, atoms = frame
                
                # 跳过空帧
                if n_atoms_orig == 0:
                    continue
                
                # 计算几何中心
                sum_x = sum(atom[1] for atom in atoms)
                sum_y = sum(atom[2] for atom in atoms)
                sum_z = sum(atom[3] for atom in atoms)
                xc = sum_x / n_atoms_orig
                yc = sum_y / n_atoms_orig
                zc = sum_z / n_atoms_orig
                
                # 计算距离并排序
                distances = [
                    ( (atom[1]-xc)**2 + (atom[2]-yc)**2 + (atom[3]-zc)**2, idx )
                    for idx, atom in enumerate(atoms)
                ]
                distances.sort()
                
                # 确定实际要删除的数量（不能超过原子总数）
                actual_remove = min(count, n_atoms_orig)
                
                # 获取要删除的原子索引
                remove_indices = { idx for _, idx in distances[:actual_remove] }
                
                # 保留未删除的原子
                remaining_atoms = [atom for idx, atom in enumerate(atoms) if idx not in remove_indices]
                
                # 写入处理后的帧
                outfile.write(f"{len(remaining_atoms)}\n")
                outfile.write(f"{comment}\n")
                for atom in remaining_atoms:
                    outfile.write(f"{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n")

if __name__ == "__main__":
    input_filename = "total.xyz"  # 输入文件路径
    remove_counts = [1, 2, 4, 6, 8, 10]  # 需要删除的原子数量列表
    
    # 检查输入文件是否存在
    if not os.path.exists(input_filename):
        print(f"错误：输入文件 {input_filename} 不存在！")
        sys.exit(1)
    
    # 执行处理
    process_frames(input_filename, remove_counts)
    print(f"处理完成！生成以下文件：")
    for count in remove_counts:
        print(f"  output_{count}.xyz")
