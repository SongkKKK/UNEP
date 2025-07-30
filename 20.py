def filter_exyz_frames(input_exyz, output_interesting, output_uninteresting):
    # 定义感兴趣的元素集合
    interesting_elements = {'Ag', 'Al', 'Au', 'Cr', 'Cu', 'Mg', 'Mo', 'Ni', 'Pb', 'Pd', 'Pt', 'Ta', 'Ti', 'V', 'W', 'Zr', 'Hf', 'Nb', 'Re', 'Os'}    
    
    with open(input_exyz, 'r') as infile, \
         open(output_interesting, 'w') as out_interesting, \
         open(output_uninteresting, 'w') as out_uninteresting:
        
        while True:
            # 读取并解析帧的头部（原子总数和注释行）
            atom_count_line = infile.readline()
            if not atom_count_line:  # 检查文件末尾
                break
                
            comment_line = infile.readline()
            
            # 尝试转换原子总数
            try:
                atom_count = int(atom_count_line.strip())
            except ValueError:
                continue
            
            # 读取帧中的原子数据
            atoms = [infile.readline() for _ in range(atom_count)]
            
            # 检查帧中的原子是否全部属于感兴趣的元素
            if all(atom.split()[0] in interesting_elements for atom in atoms):
                # 符合条件 - 写入interesting文件
                out_interesting.write(f"{atom_count}\n")
                out_interesting.write(comment_line)
                out_interesting.writelines(atoms)
            else:
                # 不符合条件 - 写入uninteresting文件
                out_uninteresting.write(f"{atom_count}\n")
                out_uninteresting.write(comment_line)
                out_uninteresting.writelines(atoms)

# 脚本使用示例
input_exyz = 'DEEP2XYZ.xyz'           # 输入文件路径
output_interesting = '20.xyz'      # 符合条件的输出文件
output_uninteresting = 'uninteresting.xyz'  # 不符合条件的输出文件

filter_exyz_frames(input_exyz, output_interesting, output_uninteresting)
