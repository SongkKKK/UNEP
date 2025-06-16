from ase.build import bulk
from ase.io import write
import random
import os

# 定义输出目录
output_dir = 'hcp'
os.makedirs(output_dir, exist_ok=True)

# 定义所有元素及其组合
elements = ["Ag", "Al", "Au", "Cr", "Cu", "Mg", "Mo", "Ni", "Pb", "Pd", "Pt", "Ta", "Ti", "V", "W", "Zr", "Hf", "Nb", "Os", "Re"]
# 通过组合生成元素对
element_pairs = [(elements[i], elements[j]) for i in range(len(elements)) for j in range(i+1, len(elements))]

# 浓度
concentrations = [0.02, 0.05, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
supercell_size = (3, 3, 2)

for elem1, elem2 in element_pairs:
    for conc in concentrations:
        # 创建BCC基础结构
        #base_atoms = bulk(elem1, 'bcc', a=3.5, cubic=True).repeat(supercell_size)
        base_atoms = bulk(elem1, 'hcp', a=4, c=5.5).repeat(supercell_size)
        num_atoms = len(base_atoms)
        
        # 计算第二种元素的原子数
        num_elem2_atoms = int(conc * num_atoms)
        
        # 随机选择原子进行替换
        indices = random.sample(range(num_atoms), num_elem2_atoms)
        symbols = base_atoms.get_chemical_symbols()
        
        for index in indices:
            symbols[index] = elem2
        
        # 更新原子的符号
        base_atoms.set_chemical_symbols(symbols)
        
        # 文件命名
        conc_str = str(conc).replace('.', 'p')
        output_file = os.path.join(output_dir, f'{elem1}_{elem2}_bcc_{conc_str}.xyz')
        
        # 输出结构到文件
        write(output_file, base_atoms)

print('二元BCC结构已生成并保存到对应目录下。')
