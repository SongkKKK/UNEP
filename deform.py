from wizard.atoms import SymbolInfo, Morph
from wizard.frames import MultiMol
from wizard.io import relax, read_xyz

frames = read_xyz('init.xyz')
Temp_Struts = []
temperatures = [300, 500]
for i, atoms in enumerate(frames):
    for temperature in temperatures:
        dirname = f'deform/{i}/{temperature}K/relax'
        strain_rate = 1e9 #可修改
        dt = 1e-15 #fs
        length = atoms.cell[2, 2]
        strain = strain_rate * dt * length
        run_in=['potential nep.txt',
                f'velocity {temperature}',   
                'time_step 1',
                f'ensemble npt_ber {temperature} {temperature} 100 0 0 0 100 100 100 1000',
                f'deform {strain} 0 0 1',
                #f'mc canonical 100 100 {temperature} {temperature}',
                'dump_thermo 1000',
                'dump_exyz 10000',
                'run 500000']
        Morph(atoms).gpumd(dirname=dirname, run_in=run_in)
        Temp_Struts += read_xyz(dirname + '/dump.xyz')

MultiMol(Temp_Struts).dump('deform.xyz')
