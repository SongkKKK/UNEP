import dpdata
perturbed_system = dpdata.System('./POSCAR').perturb(pert_num=60, 
    cell_pert_fraction=0.1, 
    atom_pert_distance=0.2, 
    atom_pert_style='normal')
print(perturbed_system.data)
for i in range(60):
	perturbed_system.to_vasp_poscar('POSCARs/POSCAR%s'%i,frame_idx=i)
