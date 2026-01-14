import numpy as np
import random
from itertools import combinations
import os
from tqdm import tqdm

# 参数设置
elements = ["Ag", "Al", "Au", "Cr", "Cu", "Mg", "Mo", "Ni", "Pb", "Pd"]
lattice_constants = {
    "bcc": 3.8,  # BCC晶格常数 (Å)
    "fcc": 3.8,  # FCC晶格常数 (Å)
    "hcp": {"a": 3.2, "c": 5.2}  # HCP晶格常数 (Å)
}

# 不同结构的超胞尺寸和单位晶胞原子数
structures_config = {
    "bcc": {
        "supercell_size": [3, 3, 3],
        "atoms_per_unit": 2,  # BCC单位晶胞原子数
        "description": "BCC (Body-Centered Cubic)"
    },
    "fcc": {
        "supercell_size": [2, 2, 2],
        "atoms_per_unit": 4,  # FCC单位晶胞原子数
        "description": "FCC (Face-Centered Cubic)"
    },
    "hcp": {
        "supercell_size": [3, 3, 2],
        "atoms_per_unit": 2,  # HCP单位晶胞原子数
        "description": "HCP (Hexagonal Close-Packed)"
    }
}

concentration_step = 0.05  # 浓度变化间隔
min_concentration = 0.05  # 最小浓度
sampling_step = 2  # 每隔2个浓度点采样一次

# 创建输出目录
output_dir = "binary_alloys_multistructure"
os.makedirs(output_dir, exist_ok=True)

# 生成BCC超胞坐标（分数坐标）
def generate_bcc_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                # BCC晶胞内的两个原子位置 (分数坐标)
                positions.append([i, j, k])            # 角点
                positions.append([i + 0.5, j + 0.5, k + 0.5])  # 体心
    return np.array(positions)

# 生成FCC超胞坐标（分数坐标）
def generate_fcc_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                # FCC晶胞内的四个原子位置 (分数坐标)
                positions.append([i, j, k])               # 角点 (0,0,0)
                positions.append([i + 0.5, j + 0.5, k])    # 面心 (1/2,1/2,0)
                positions.append([i + 0.5, j, k + 0.5])    # 面心 (1/2,0,1/2)
                positions.append([i, j + 0.5, k + 0.5])    # 面心 (0,1/2,1/2)
    return np.array(positions)

# 生成HCP超胞坐标（分数坐标）
def generate_hcp_positions(size_x, size_y, size_z):
    positions = []
    a = lattice_constants["hcp"]["a"]
    c = lattice_constants["hcp"]["c"]
    
    # HCP晶胞参数
    a1 = np.array([a, 0, 0])
    a2 = np.array([-a/2, a*np.sqrt(3)/2, 0])
    a3 = np.array([0, 0, c])
    
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                # 六方晶胞内的两个原子位置 (分数坐标)
                # 位置1: (0, 0, 0)
                positions.append([i, j, k])
                # 位置2: (1/3, 2/3, 1/2)
                positions.append([i + 1/3, j + 2/3, k + 0.5])
    
    return np.array(positions)

# 将分数坐标转换为直角坐标
def fractional_to_cartesian(frac_pos, structure_type):
    if structure_type == "hcp":
        # HCP需要特殊处理，因为它不是立方晶系
        a = lattice_constants["hcp"]["a"]
        c = lattice_constants["hcp"]["c"]
        
        # HCP的基矢量
        a1 = np.array([a, 0, 0])
        a2 = np.array([-a/2, a*np.sqrt(3)/2, 0])
        a3 = np.array([0, 0, c])
        
        # 转换为直角坐标
        cart_pos = np.zeros((len(frac_pos), 3))
        for idx, pos in enumerate(frac_pos):
            cart_pos[idx] = pos[0]*a1 + pos[1]*a2 + pos[2]*a3
        
        # 超胞尺寸
        cell_size_x = a * structures_config["hcp"]["supercell_size"][0]
        cell_size_y = a * np.sqrt(3)/2 * structures_config["hcp"]["supercell_size"][1]
        cell_size_z = c * structures_config["hcp"]["supercell_size"][2]
        
        # 应用周期性边界条件
        cart_pos[:, 0] = np.mod(cart_pos[:, 0], cell_size_x)
        cart_pos[:, 1] = np.mod(cart_pos[:, 1], cell_size_y)
        cart_pos[:, 2] = np.mod(cart_pos[:, 2], cell_size_z)
        
        return cart_pos, (cell_size_x, cell_size_y, cell_size_z)
    else:
        # BCC和FCC的处理方式相同（立方晶系）
        lattice_const = lattice_constants[structure_type]
        
        # 超胞尺寸
        supercell_size = structures_config[structure_type]["supercell_size"]
        cell_size_x = lattice_const * supercell_size[0]
        cell_size_y = lattice_const * supercell_size[1]
        cell_size_z = lattice_const * supercell_size[2]
        
        # 转换为直角坐标
        cart_pos = frac_pos * lattice_const
        
        # 应用周期性边界条件
        cart_pos[:, 0] = np.mod(cart_pos[:, 0], cell_size_x)
        cart_pos[:, 1] = np.mod(cart_pos[:, 1], cell_size_y)
        cart_pos[:, 2] = np.mod(cart_pos[:, 2], cell_size_z)
        
        return cart_pos, (cell_size_x, cell_size_y, cell_size_z)

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
def generate_binary_alloy(elements, concentration, positions, structure_type):
    # 确保浓度有效
    if concentration < min_concentration or concentration > (1.0 - min_concentration):
        return None, None
    
    # 计算两种元素的浓度
    conc1 = concentration
    conc2 = 1.0 - concentration
    
    # 获取总原子数
    config = structures_config[structure_type]
    total_atoms = config["supercell_size"][0] * config["supercell_size"][1] * config["supercell_size"][2] * config["atoms_per_unit"]
    
    # 计算原子数量
    n1 = max(1, round(total_atoms * conc1))
    n2 = total_atoms - n1
    
    # 确保原子数量有效
    if n1 < 1 or n2 < 1:
        return None, None
    
    # 创建原子类型列表
    atom_types = [elements[0]] * n1 + [elements[1]] * n2
    random.shuffle(atom_types)  # 随机分布
    
    # 转换为直角坐标并获取晶胞尺寸
    cart_positions, cell_sizes = fractional_to_cartesian(positions, structure_type)
    
    # 实际浓度
    actual_conc1 = n1 / total_atoms
    actual_conc2 = n2 / total_atoms
    
    # 生成XYZ格式内容
    xyz_content = f"{total_atoms}\n"
    
    if structure_type == "hcp":
        # HCP需要特殊的晶格矢量表示
        a = lattice_constants["hcp"]["a"]
        c = lattice_constants["hcp"]["c"]
        x_size, y_size, z_size = cell_sizes
        xyz_content += f"Lattice=\"{x_size} 0 0 {y_size/2} {y_size*np.sqrt(3)/2} 0 0 0 {z_size}\" "
    else:
        # BCC和FCC使用简单的立方晶格表示
        x_size, y_size, z_size = cell_sizes
        xyz_content += f"Lattice=\"{x_size} 0 0 0 {y_size} 0 0 0 {z_size}\" "
    
    xyz_content += f"Properties=species:S:1:pos:R:3 "
    xyz_content += f"Structure={structure_type.upper()} "
    xyz_content += f"Binary={elements[0]}-{elements[1]} "
    xyz_content += f"Target_Conc={conc1:.2f}-{conc2:.2f} "
    xyz_content += f"Actual_Conc={actual_conc1:.4f}-{actual_conc2:.4f} "
    xyz_content += f"Supercell={config['supercell_size'][0]}x{config['supercell_size'][1]}x{config['supercell_size'][2]}\n"
    
    for atom, pos in zip(atom_types, cart_positions):
        xyz_content += f"{atom} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n"
    
    return xyz_content, total_atoms

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
    
    # 计算总结构数
    total_structures_all = 0
    structures_info = {}
    
    # 为每种结构类型生成位置和计算总原子数
    structure_positions = {}
    structure_totals = {}
    
    for structure_type, config in structures_config.items():
        supercell_size = config["supercell_size"]
        
        # 生成原子位置
        if structure_type == "bcc":
            positions = generate_bcc_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        elif structure_type == "fcc":
            positions = generate_fcc_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        elif structure_type == "hcp":
            positions = generate_hcp_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        
        structure_positions[structure_type] = positions
        
        # 计算总原子数
        total_atoms = supercell_size[0] * supercell_size[1] * supercell_size[2] * config["atoms_per_unit"]
        structure_totals[structure_type] = total_atoms
        
        # 计算该结构的总结构数
        structures_per_type = len(element_combinations) * len(concentrations)
        structures_info[structure_type] = structures_per_type
        total_structures_all += structures_per_type
        
        print(f"\n{config['description']}:")
        print(f"  Supercell: {supercell_size[0]}x{supercell_size[1]}x{supercell_size[2]}")
        print(f"  Atoms per cell: {config['atoms_per_unit']}")
        print(f"  Total atoms: {total_atoms}")
        print(f"  Estimated structures: {structures_per_type}")
    
    print(f"\nEstimated total structures (all types): {total_structures_all}")
    
    # 创建总进度条
    pbar = tqdm(total=total_structures_all, desc="Generating all structures")
    
    total_structures_created = 0
    
    # 为每种结构类型、每个元素组合和每个浓度生成结构
    for structure_type, config in structures_config.items():
        struct_output_dir = os.path.join(output_dir, structure_type)
        os.makedirs(struct_output_dir, exist_ok=True)
        
        positions = structure_positions[structure_type]
        
        for pair in element_combinations:
            pair_dir = os.path.join(struct_output_dir, f"{pair[0]}_{pair[1]}")
            os.makedirs(pair_dir, exist_ok=True)
            
            for conc in concentrations:
                xyz_data, atoms_count = generate_binary_alloy(
                    pair, conc, positions, structure_type
                )
                
                if xyz_data:
                    # 创建文件名
                    conc_str = f"{conc:.2f}".replace('.', 'p')
                    filename = f"{pair[0]}_{pair[1]}_{structure_type}_c_{conc_str}.xyz"
                    filepath = os.path.join(pair_dir, filename)
                    
                    # 保存到文件
                    with open(filepath, 'w') as f:
                        f.write(xyz_data)
                    
                    total_structures_created += 1
                
                # 更新进度条
                pbar.update(1)
    
    pbar.close()
    
    print("\nAll structures generated successfully!")
    print(f"Total structures created: {total_structures_created}")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    
    print("\nStructure details:")
    for structure_type, config in structures_config.items():
        supercell = config["supercell_size"]
        total_atoms = structure_totals[structure_type]
        print(f"  {config['description']}: {supercell[0]}x{supercell[1]}x{supercell[2]} ({total_atoms} atoms)")
