import numpy as np

def write_exyz_file(file_name, element, positions, lattice_vectors):
    with open(file_name, 'w') as f:
        f.write(f"{len(positions)}\n")
        lattice_str = " ".join(f"{v:.6f}" for vector in lattice_vectors for v in vector)
        f.write(f'pbc="T T T" Lattice="{lattice_str}" Properties=species:S:1:pos:R:3\n')
        for position in positions:
            f.write(f"{element} {position[0]:.6f} {position[1]:.6f} {position[2]:.6f}\n")

def generate_supercell(lattice_vectors, basis, repetitions):
    positions = []
    for x in range(repetitions[0]):
        for y in range(repetitions[1]):
            for z in range(repetitions[2]):
                translation = x * lattice_vectors[0] + y * lattice_vectors[1] + z * lattice_vectors[2]
                for atom in basis:
                    positions.append(atom + translation)
    expanded_vectors = np.array(lattice_vectors) * np.array(repetitions).reshape(-1, 1)
    return np.array(positions), expanded_vectors

def generate_bcc(a, element, repetitions):
    lattice_vectors = np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
    basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) * a
    positions, expanded_vectors = generate_supercell(lattice_vectors, basis, repetitions)
    write_exyz_file(f"{element}_bcc.exyz", element, positions, expanded_vectors)

def generate_fcc(a, element, repetitions):
    lattice_vectors = np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
    basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]]) * a
    positions, expanded_vectors = generate_supercell(lattice_vectors, basis, repetitions)
    write_exyz_file(f"{element}_fcc.exyz", element, positions, expanded_vectors)

def generate_hcp(a, c, element, repetitions):
    lattice_vectors = np.array([
        [a, 0, 0],
        [-a / 2, a * (3**0.5) / 2, 0],
        [0, 0, c]
    ])
    basis = np.array([[0.0, 0.0, 0.0], [1/3, 2/3, 0.5]]) * a
    basis[:, 2] *= c / a
    positions, expanded_vectors = generate_supercell(lattice_vectors, basis, repetitions)
    write_exyz_file(f"{element}_hcp.exyz", element, positions, expanded_vectors)

elements = ["Al", "Mg", "Zn", "Zr", "Mo"]
a_bcc = a_fcc = 4.0
a_hcp = 4.0
c_hcp = 5.0

repetitions_bcc = (3, 3, 3)
repetitions_fcc = (2, 2, 2)
repetitions_hcp = (4, 4, 3)

for element in elements:
    generate_bcc(a_bcc, element, repetitions_bcc)
    generate_fcc(a_fcc, element, repetitions_fcc)
    generate_hcp(a_hcp, c_hcp, element, repetitions_hcp)
