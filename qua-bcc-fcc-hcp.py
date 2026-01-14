import numpy as np
import random
from itertools import combinations
import os
from tqdm import tqdm

# 参数设置
elements = ["Ag", "Al", "Au", "Cr", "Cu", "Mg", "Mo", "Ni", "Pb", "Pd", 
            "Pt", "Ta", "Ti", "V", "W", "Zr", "Hf", "Nb", "Os", "Re"]

# 不同结构的晶格常数和配置
structures_config = {
    "bcc": {
        "supercell_size": [3, 3, 3],
        "atoms_per_unit": 2,  # BCC单位晶胞原子数
        "lattice_constant": 3.6,  # BCC晶格常数 (Å)
        "description": "BCC (Body-Centered Cubic)",
        "lattice_type": "cubic"
    },
    "fcc": {
        "supercell_size": [2, 2, 2],
        "atoms_per_unit": 4,  # FCC单位晶胞原子数
        "lattice_constant": 3.6,  # FCC晶格常数 (Å)
        "description": "FCC (Face-Centered Cubic)",
        "lattice_type": "cubic"
    },
    "hcp": {
        "supercell_size": [3, 3, 2],
        "atoms_per_unit": 2,  # HCP单位晶胞原子数
        "lattice_constant": {"a": 3.2, "c": 5.2},  # HCP晶格常数 (Å)
        "description": "HCP (Hexagonal Close-Packed)",
        "lattice_type": "hexagonal"
    }
}

concentration_step = 0.2  # 浓度变化间隔
min_concentration = 0.2   # 最小浓度
sampling_step = 1         # 采样步长

# 创建输出目录
output_dir = "quaternary_alloys_multistructure"
os.makedirs(output_dir, exist_ok=True)

# 生成BCC超胞坐标（分数坐标）
def generate_bcc_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                positions.append([i, j, k])            # 角点
                positions.append([i + 0.5, j + 0.5, k + 0.5])  # 体心
    return np.array(positions)

# 生成FCC超胞坐标（分数坐标）
def generate_fcc_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                positions.append([i, j, k])               # 角点
                positions.append([i, j + 0.5, k + 0.5])    # 面心1
                positions.append([i + 0.5, j, k + 0.5])    # 面心2
                positions.append([i + 0.5, j + 0.5, k])    # 面心3
    return np.array(positions)

# 生成HCP超胞坐标（分数坐标）- 六方晶系
def generate_hcp_positions(size_x, size_y, size_z):
    positions = []
    for i in range(size_x):
        for j in range(size_y):
            for k in range(size_z):
                # HCP晶胞内的两个原子位置
                positions.append([i, j, k])                # 位置1
                positions.append([i + 1/3, j + 2/3, k + 0.5])  # 位置2
    return np.array(positions)

# 将分数坐标转换为直角坐标
def fractional_to_cartesian(frac_pos, structure_type, config):
    if structure_type == "hcp":
        # HCP六方晶系转换
        a = config["lattice_constant"]["a"]
        c = config["lattice_constant"]["c"]
        supercell_size = config["supercell_size"]
        
        # 六方晶系基矢
        a1 = np.array([a, 0, 0])
        a2 = np.array([-a/2, a*np.sqrt(3)/2, 0])
        a3 = np.array([0, 0, c])
        
        # 计算超胞尺寸
        cell_size_x = a * supercell_size[0]
        cell_size_y = a * np.sqrt(3)/2 * supercell_size[1]
        cell_size_z = c * supercell_size[2]
        
        # 转换为直角坐标
        cart_pos = np.zeros((len(frac_pos), 3))
        for idx, pos in enumerate(frac_pos):
            cart_pos[idx] = pos[0]*a1 + pos[1]*a2 + pos[2]*a3
        
        # 应用周期性边界条件
        cart_pos[:, 0] = np.mod(cart_pos[:, 0], cell_size_x)
        cart_pos[:, 1] = np.mod(cart_pos[:, 1], cell_size_y)
        cart_pos[:, 2] = np.mod(cart_pos[:, 2], cell_size_z)
        
        return cart_pos, (cell_size_x, cell_size_y, cell_size_z)
    else:
        # BCC和FCC立方晶系转换
        lattice_const = config["lattice_constant"]
        supercell_size = config["supercell_size"]
        
        # 计算超胞尺寸
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
    ideal_counts = [c * total_atoms for c in concentrations]
    counts = [max(1, round(count)) for count in ideal_counts]
    total = sum(counts)
    
    if total != total_atoms:
        fractional_parts = [count - int(count) for count in ideal_counts]
        indices = np.argsort(fractional_parts)[::-1]
        
        while total != total_atoms:
            if total > total_atoms:
                idx = indices[-1]
                if counts[idx] > 1:
                    counts[idx] -= 1
                    total -= 1
                indices = indices[:-1]
                if len(indices) == 0:
                    indices = np.argsort(fractional_parts)[::-1]
            else:
                idx = indices[0]
                counts[idx] += 1
                total += 1
                indices = indices[1:]
                if len(indices) == 0:
                    indices = np.argsort(fractional_parts)[::-1]
    
    return counts

# 生成四元合金结构
def generate_quaternary_alloy(elements, concentrations, positions, structure_type):
    # 检查浓度有效性
    if abs(sum(concentrations) - 1.0) > 1e-6:
        return None
    if any(c < min_concentration for c in concentrations):
        return None
    
    # 获取结构配置
    config = structures_config[structure_type]
    
    # 计算总原子数
    total_atoms = (config["supercell_size"][0] * 
                   config["supercell_size"][1] * 
                   config["supercell_size"][2] * 
                   config["atoms_per_unit"])
    
    # 计算原子数量
    atom_counts = calculate_atom_counts(concentrations, total_atoms)
    
    # 创建原子类型列表
    atom_types = []
    for i, elem in enumerate(elements):
        atom_types.extend([elem] * atom_counts[i])
    
    # 随机分布
    random.shuffle(atom_types)
    
    # 转换为直角坐标
    cart_positions, cell_sizes = fractional_to_cartesian(positions, structure_type, config)
    
    # 计算实际浓度
    actual_concentrations = [count / total_atoms for count in atom_counts]
    
    # 生成XYZ格式内容
    xyz_content = f"{total_atoms}\n"
    
    # 添加晶格信息
    if structure_type == "hcp":
        # HCP特殊格式
        a = config["lattice_constant"]["a"]
        c = config["lattice_constant"]["c"]
        xyz_content += f"Lattice=\"{cell_sizes[0]:.4f} 0 0 0 {cell_sizes[1]:.4f} 0 0 0 {cell_sizes[2]:.4f}\" "
        xyz_content += f"HCP_params=\"{a:.4f}_{c:.4f}\" "
        xyz_content += f"c/a_ratio=\"{c/a:.4f}\" "
    else:
        # BCC和FCC格式
        xyz_content += f"Lattice=\"{cell_sizes[0]:.4f} 0 0 0 {cell_sizes[1]:.4f} 0 0 0 {cell_sizes[2]:.4f}\" "
        xyz_content += f"Lattice_constant=\"{config['lattice_constant']:.4f}\" "
    
    xyz_content += f"Properties=species:S:1:pos:R:3 "
    xyz_content += f"Structure={structure_type.upper()} "
    xyz_content += f"Quaternary={'-'.join(elements)} "
    xyz_content += f"Target_Conc={'-'.join([f'{c:.2f}' for c in concentrations])} "
    xyz_content += f"Actual_Conc={'-'.join([f'{ac:.4f}' for ac in actual_concentrations])} "
    xyz_content += f"Supercell={config['supercell_size'][0]}x{config['supercell_size'][1]}x{config['supercell_size'][2]} "
    xyz_content += f"Atoms_per_cell={config['atoms_per_unit']}\n"
    
    for atom, pos in zip(atom_types, cart_positions):
        xyz_content += f"{atom} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n"
    
    return xyz_content

# 生成稀疏采样的浓度组合
def generate_sparse_concentration_combinations():
    concentrations = []
    
    # 生成所有可能的浓度值
    values = [round(i * concentration_step, 2) for i in 
              range(int(min_concentration / concentration_step), 
                    int((1.0 - 3*min_concentration) / concentration_step) + 1)]
    
    # 稀疏采样：每隔sampling_step个点取一个
    sparse_values = values[::sampling_step]
    
    # 生成四元浓度组合（稀疏采样）
    total_combinations = 0
    for c1 in sparse_values:
        for c2 in sparse_values:
            for c3 in sparse_values:
                c4 = 1.0 - c1 - c2 - c3
                if c4 >= min_concentration:
                    concentrations.append((c1, c2, c3, c4))
                    total_combinations += 1
    
    return concentrations

# 主程序
if __name__ == "__main__":
    # 生成所有可能的四元组合
    element_combinations = list(combinations(elements, 4))
    print(f"Total element combinations: {len(element_combinations)}")
    
    # 生成稀疏采样的浓度组合
    concentration_combinations = generate_sparse_concentration_combinations()
    print(f"Total sparse concentration combinations: {len(concentration_combinations)}")
    
    # 显示结构配置信息
    print("\n" + "="*60)
    print("STRUCTURE CONFIGURATIONS")
    print("="*60)
    
    total_structures_all = 0
    structure_details = {}
    
    for structure_type, config in structures_config.items():
        supercell = config["supercell_size"]
        total_atoms = supercell[0] * supercell[1] * supercell[2] * config["atoms_per_unit"]
        
        # 计算该结构的配置信息
        structures_per_type = len(element_combinations) * len(concentration_combinations)
        total_structures_all += structures_per_type
        
        structure_details[structure_type] = {
            "total_atoms": total_atoms,
            "structures_per_type": structures_per_type
        }
        
        print(f"\n{config['description']}:")
        print(f"  Supercell: {supercell[0]}x{supercell[1]}x{supercell[2]}")
        print(f"  Atoms per unit cell: {config['atoms_per_unit']}")
        print(f"  Total atoms in supercell: {total_atoms}")
        print(f"  Lattice type: {config['lattice_type']}")
        
        if structure_type == "hcp":
            a = config["lattice_constant"]["a"]
            c = config["lattice_constant"]["c"]
            print(f"  Lattice parameters: a={a:.4f} Å, c={c:.4f} Å")
            print(f"  c/a ratio: {c/a:.4f}")
        else:
            print(f"  Lattice constant: {config['lattice_constant']:.4f} Å")
        
        print(f"  Estimated structures: {structures_per_type:,}")
    
    print(f"\nEstimated total structures (all types): {total_structures_all:,}")
    print("="*60)
    
    # 预先生成每种结构的位置
    structure_positions = {}
    for structure_type, config in structures_config.items():
        supercell_size = config["supercell_size"]
        
        if structure_type == "bcc":
            positions = generate_bcc_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        elif structure_type == "fcc":
            positions = generate_fcc_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        elif structure_type == "hcp":
            positions = generate_hcp_positions(supercell_size[0], supercell_size[1], supercell_size[2])
        
        structure_positions[structure_type] = positions
    
    # 创建总进度条
    pbar = tqdm(total=total_structures_all, desc="Generating all quaternary structures")
    total_structures_created = 0
    
    # 为每种结构类型生成结构
    for structure_type, config in structures_config.items():
        struct_output_dir = os.path.join(output_dir, structure_type)
        os.makedirs(struct_output_dir, exist_ok=True)
        
        positions = structure_positions[structure_type]
        
        print(f"\nGenerating {config['description']} structures...")
        
        for quartet in element_combinations:
            quartet_dir = os.path.join(struct_output_dir, f"{'_'.join(quartet)}")
            os.makedirs(quartet_dir, exist_ok=True)
            
            for conc in concentration_combinations:
                concentrations = list(conc)
                xyz_data = generate_quaternary_alloy(
                    quartet, concentrations, positions, structure_type
                )
                
                if xyz_data:
                    # 创建文件名
                    conc_str = '_'.join([f"{c:.2f}" for c in concentrations]).replace('.', 'p')
                    filename = f"{'_'.join(quartet)}_{structure_type}_c_{conc_str}.xyz"
                    filepath = os.path.join(quartet_dir, filename)
                    
                    # 保存到文件
                    with open(filepath, 'w') as f:
                        f.write(xyz_data)
                    
                    total_structures_created += 1
                
                # 更新进度条
                pbar.update(1)
    
    pbar.close()
    
    print("\n" + "="*60)
    print("GENERATION COMPLETE")
    print("="*60)
    
    print(f"\nTotal structures created: {total_structures_created:,}")
    print(f"Output directory: {os.path.abspath(output_dir)}")
    
    print("\nStructure summary:")
    for structure_type, config in structures_config.items():
        details = structure_details[structure_type]
        supercell = config["supercell_size"]
        print(f"  {config['description']}:")
        print(f"    Supercell: {supercell[0]}x{supercell[1]}x{supercell[2]}")
        print(f"    Total atoms: {details['total_atoms']}")
        print(f"    Estimated structures: {details['structures_per_type']:,}")
    
    print("\nStatistics:")
    print(f"  Number of elements: {len(elements)}")
    print(f"  Number of element quartets: {len(element_combinations):,}")
    print(f"  Number of concentration combinations: {len(concentration_combinations):,}")
    print(f"  Concentration step: {concentration_step}")
    print(f"  Minimum concentration: {min_concentration}")
    print(f"  Sampling step: {sampling_step}")
    
    print("\nFile organization:")
    print(f"  {output_dir}/")
    for structure_type in structures_config.keys():
        print(f"    {structure_type}/")
        print(f"      element1_element2_element3_element4/")
        print(f"        *_{structure_type}_c_*.xyz files")
    
    print("\nNote: This script generates quaternary alloy structures with")
    print(f"      {concentration_step} concentration steps and {min_concentration} minimum concentration.")
    print("      Each structure is randomly shuffled to ensure proper mixing.")