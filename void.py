import sys

def remove_atoms_from_center(input_file, output_file, remove_positions_1based):
    # 将1-based的删除位置转换为0-based索引
    remove_positions = [pos - 1 for pos in remove_positions_1based]
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            # 读取原子数目行
            line = infile.readline()
            if not line:
                break  # 文件结束
            n_atoms = int(line.strip())
            
            # 读取注释行
            comment = infile.readline().strip()
            
            # 读取原子数据
            atoms = []
            for _ in range(n_atoms):
                parts = infile.readline().split()
                if len(parts) < 4:
                    print("错误：原子行格式不正确。")
                    sys.exit(1)
                elem, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
                atoms.append((elem, x, y, z))
            
            # 计算几何中心
            if n_atoms == 0:
                continue  # 跳过空帧
            sum_x = sum(atom[1] for atom in atoms)
            sum_y = sum(atom[2] for atom in atoms)
            sum_z = sum(atom[3] for atom in atoms)
            xc = sum_x / n_atoms
            yc = sum_y / n_atoms
            zc = sum_z / n_atoms
            
            # 计算每个原子到中心的距离平方并排序
            distances = []
            for idx, atom in enumerate(atoms):
                dx = atom[1] - xc
                dy = atom[2] - yc
                dz = atom[3] - zc
                dist_sq = dx**2 + dy**2 + dz**2
                distances.append((dist_sq, idx))
            
            # 按距离从小到大排序
            distances.sort()
            
            # 获取要删除的原子索引（处理超出范围的情况）
            sorted_indices = [idx for (_, idx) in distances]
            remove_indices = []
            for pos in remove_positions:
                if pos < len(sorted_indices):
                    remove_indices.append(sorted_indices[pos])
            
            # 移除重复的索引（如果有）
            remove_indices = list(set(remove_indices))
            
            # 保留未被删除的原子
            remaining_atoms = [atom for idx, atom in enumerate(atoms) if idx not in remove_indices]
            
            # 写入处理后的帧
            outfile.write(f"{len(remaining_atoms)}\n")
            outfile.write(f"{comment}\n")
            for atom in remaining_atoms:
                outfile.write(f"{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n")

# 使用示例
if __name__ == "__main__":
    input_filename = "total.xyz"    # 输入文件名
    output_filename = "output.xyz"  # 输出文件名
    # 要删除的原子位置（1-based，即第1、2、4、6、8、10近的原子）
    remove_positions = [1, 2, 4, 6, 8, 10]
    
    remove_atoms_from_center(input_filename, output_filename, remove_positions)
    print("处理完成！")
