from ase.io import read, write
from ase import Atoms
from ase.build import bulk

def get_oct_sites_bcc(cell):
    # 返回BCC的八面体间隙位置
    return [
        (0.5, 0.0, 0.5), # 顶点到边中心
        (0.5, 0.5, 0.0), # 面中心
        (0.0, 0.5, 0.5)  # 面中心
    ]

def get_oct_sites_fcc(cell):
    # 返回FCC的八面体间隙位置
    return [
        (0.5, 0.5, 0.5), # 体心位置，仅简单示例
        (0.5, 0.0, 0.0), # 边中心
        (0.0, 0.5, 0.0), # 边中心
        (0.0, 0.0, 0.5)  # 边中心
    ]

def get_oct_sites_hcp(cell):
    # 返回HCP的八面体间隙位置（简化处理）
    return [
        (1/3, 1/3, 0.5), # 沿c轴中心
        (2/3, 2/3, 0.5)  # 对称位置
    ]

def insert_atoms(atoms, num_atoms, element, get_sites_fn):
    cell = atoms.get_cell()
    sites = get_sites_fn(cell)
    
    if len(sites) < num_atoms:
        raise ValueError("Not enough octahedral sites available for insertion!")
    
    new_atoms = Atoms(
        [element] * num_atoms,
        positions=[cell.dot(s) for s in sites[:num_atoms]],
        cell=cell,
        pbc=True
    )
    return atoms + new_atoms

def process_xyz(input_filename, output_filename):
    frames = read(input_filename, index=':')  # Read all frames

    doped_structures = []
    for atoms in frames:
        # 简单检查晶格类型，假设原子数不变，简化判断方式
        if len(atoms) % 2 == 0 and len(atoms) % 6 != 0:
            get_sites_fn = get_oct_sites_bcc  # For BCC
        elif len(atoms) % 4 == 0:
            get_sites_fn = get_oct_sites_fcc  # For FCC
        else:
            get_sites_fn = get_oct_sites_hcp  # For HCP
        
        for n in range(1, 5):
            for element in ['Al', 'Mg', 'Hf']:
                try:
                    modified_atoms = insert_atoms(atoms.copy(), n, element, get_sites_fn)
                    doped_structures.append(modified_atoms)
                except ValueError as e:
                    print(f"Skipping: {str(e)}")
    
    write(output_filename, doped_structures)

# Use your actual filenames
input_xyz = 'input.xyz'
output_xyz = 'doped_output.xyz'

process_xyz(input_xyz, output_xyz)
