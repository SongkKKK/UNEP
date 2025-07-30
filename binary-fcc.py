import numpy as np
import random
from itertools import combinations
import os
from tqdm import tqdm

# 参数设置
elements = ["Ag", "Al", "Au", "Cr", "Cu", "Mg", "Mo", "Ni", "Pb", "Pd"]
lattice_constant = 3.8  # 晶格常数 (Å)
supercell_size = [3, 3, 3]  # 超胞尺寸
atoms_per_unit = 4  # FCC单位晶胞原子数
total_atoms = supercell_size[0] * supercell_size[1] * supercell_size[2] * atoms_per_unit
concentration_step = 0.05  # 浓度变化间隔
min_concentration = 0.05  # 最小浓度
sampling_step = 2  # 每隔2个浓度点采样一次

# 计算超胞尺寸（Å）
cell_size_x = lattice_constant * supercell_size[0]
cell_size_y = lattice_constant * supercell_size[1]
cell_size_z = lattice_constant * supercell_size[2]

# 创建输出目录
output_dir = "fcc_binary_alloys_3x4x5"
os.makedirs(output_dir, exist_ok=True)

# 生成FCC超胞坐标（分数坐标）
def generate_fcc_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                # FCC晶胞内的四个原子位置 (分数坐标)
                positions.append([i, j, k])            # 角点
                positions.append([i, j+0.5, k+0.5])    # 面心1
                positions.append([i+0.5, j, k+0.5])    # 面心2
                positions.append([i+0.5, j+0.5, k])    # 面心3
    return np.array(positions)

# 将分数坐标转换为直角坐标
def fractional_to_cartesian(frac_pos, lattice_const):
    # 转换为直角坐标
    cart_pos = frac_pos * lattice_const
    # 应用周期性边界条件
    cart_pos[:, 0] = np.mod(cart_pos[:, 0], cell_size_x)
    cart_pos[:, 1] = np.mod(cart_pos[:, 1], cell_size_y)
    cart_pos[:, 2] = np.mod(cart_pos[:, 2], cell_size_z)
    return cart_pos

# 计算原子数量并确保有效
def calculate_atom_counts(concentrations, total_atoms):
    # 计算每种元素的理想原子数量
    ideal_counts = [c * total_atoms for c in concentrations]
    
    # 取整并确保至少有一个原子
    counts = [max(1, round(count)) for count in ideal_counts]
    
    # 调整总数
    total = sum(counts)
    if total != total_atoms:
        # 按小数部分大小排序
        fractional_parts = [count - int(count) for count in ideal_counts]
        indices = np.argsort(fractional_parts)[::-1]  # 从大到小排序
        
        # 调整原子数量
        while total != total_atoms:
            if total > total_atoms:
                # 减少小数部分最小的元素
                idx = indices[-1]
                if counts[idx] > 1:
                    counts[idx] -= 1
                    total -= 1
                # 更新索引
                indices = indices[:-1]
                if len(indices) == 0:
                    indices = np.argsort(fractional_parts)[::-1]
            else:
                # 增加小数部分最大的元素
                idx = indices[0]
                counts[idx] += 1
                total += 1
                # 更新索引
                indices = indices[1:]
                if len(indices) == 0:
                    indices = np.argsort(fractional_parts)[::-1]
    
    return counts

# 生成二元合金结构
def generate_binary_alloy(elements, concentration, positions, lattice_const):
    # 确保浓度有效
    if concentration < min_concentration or concentration > (1.0 - min_concentration):
        return None
    
    # 计算两种元素的浓度
    conc1 = concentration
    conc2 = 1.0 - concentration
    
    # 计算原子数量
    n1 = max(1, round(total_atoms * conc1))
    n2 = total_atoms - n1
    
    # 确保原子数量有效
    if n1 < 1 or n2 < 1:
        return None
    
    # 创建原子类型列表
    atom_types = [elements[0]] * n1 + [elements[1]] * n2
    random.shuffle(atom_types)  # 随机分布
    
    # 转换为直角坐标
    cart_positions = fractional_to_cartesian(positions, lattice_const)
    
    # 实际浓度
    actual_conc1 = n1 / total_atoms
    actual_conc2 = n2 / total_atoms
    
    # 生成XYZ格式内容
    xyz_content = f"{total_atoms}\n"
    xyz_content += f"Lattice=\"{cell_size_x} 0 0 0 {cell_size_y} 0 0 0 {cell_size_z}\" "
    xyz_content += f"Properties=species:S:1:pos:R:3 "
    xyz_content += f"Binary={elements[0]}-{elements[1]} "
    xyz_content += f"Target_Conc={conc1:.2f}-{conc2:.2f} "
    xyz_content += f"Actual_Conc={actual_conc1:.4f}-{actual_conc2:.4f}\n"
    
    for atom, pos in zip(atom_types, cart_positions):
        xyz_content += f"{atom} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n"
    
    return xyz_content

# 生成稀疏采样的浓度点
def generate_sparse_concentrations():
    # 生成所有可能的浓度值
    values = [round(i * concentration_step, 2) for i in 
              range(int(min_concentration / concentration_step), 
                    int((1.0 - min_concentration) / concentration_step) + 1)]
    
    # 稀疏采样：每隔sampling_step个点取一个
    sparse_values = values[::sampling_step]
    
    return sparse_values

# 主程序
if __name__ == "__main__":
    # 生成所有可能的二元组合
    element_combinations = list(combinations(elements, 2))
    print(f"Total element combinations: {len(element_combinations)}")
    
    # 生成稀疏采样的浓度点
    concentrations = generate_sparse_concentrations()
    print(f"Total sparse concentration points: {len(concentrations)}")
    estimated_total = len(element_combinations) * len(concentrations)
    print(f"Estimated total structures: {estimated_total}")
    
    # 生成FCC位置坐标
    fcc_positions = generate_fcc_positions(supercell_size[0], supercell_size[1], supercell_size[2])
    print(f"Total atoms: {total_atoms}")
    print(f"Supercell dimensions: {cell_size_x} Å × {cell_size_y} Å × {cell_size_z} Å")
    print(f"Lattice constant: {lattice_constant} Å")
    print(f"Sparse sampling: every {sampling_step} concentration points")
    
    # 创建总进度条
    total_structures = 0
    pbar_total = len(element_combinations) * len(concentrations)
    pbar = tqdm(total=pbar_total, desc="Generating structures")
    
    # 为每个组合和浓度生成结构
    for pair in element_combinations:
        pair_dir = os.path.join(output_dir, f"{pair[0]}_{pair[1]}")
        os.makedirs(pair_dir, exist_ok=True)
        
        for conc in concentrations:
            xyz_data = generate_binary_alloy(
                pair, conc, fcc_positions, lattice_constant
            )
            
            if xyz_data:
                # 创建文件名
                conc_str = f"{conc:.2f}".replace('.', 'p')
                filename = f"{pair[0]}_{pair[1]}_c_{conc_str}.xyz"
                filepath = os.path.join(pair_dir, filename)
                
                # 保存到文件
                with open(filepath, 'w') as f:
                    f.write(xyz_data)
                
                total_structures += 1
            
            # 更新进度条
            pbar.update(1)
    
    pbar.close()
    
    print("\nAll structures generated successfully!")
    print(f"Total structures created: {total_structures}")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    print(f"Supercell size: {supercell_size[0]}x{supercell_size[1]}x{supercell_size[2]} ({total_atoms} atoms)")
