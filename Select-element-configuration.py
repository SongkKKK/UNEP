from ase.io import read, write

# 设定输入和输出文件
input_filename = 'add10.xyz'  # 输入多帧XYZ文件
output_filename = 'one_structures.xyz'  # 输出包含二元元素的XYZ文件

# 读取所有帧
atoms_list = read(input_filename, index=':')

# 用于存储符合条件的帧
binary_structures = []

# 遍历所有晶胞/帧
for i, atoms in enumerate(atoms_list):
    # 获取当前帧中的元素符号，并转换为集合（去重）
    elements_in_frame = set(atoms.get_chemical_symbols())
    
    # 检查当前帧的元素种类是否正好有两种
    if len(elements_in_frame) == 1:
        print(f"Frame {i+1} is binary with elements: {elements_in_frame}")
        binary_structures.append(atoms)

# 将符合条件的帧写入输出文件
write(output_filename, binary_structures)
print(f"Binary structures have been written to {output_filename}")
