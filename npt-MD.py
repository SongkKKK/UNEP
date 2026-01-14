from wizard.atoms import SymbolInfo, Morph
from wizard.frames import MultiMol
from wizard.io import relax, read_xyz

frames = read_xyz('init.xyz')
Temp_Struts = []
temperatures = [tem1, tem2,...]
for i, atoms in enumerate(frames):
    for temperature in temperatures:
        dirname = f'npt/{i}/{temperature}K/relax'
        run_in=['potential nep.txt',
                f'velocity {temperature}',   
                'time_step 2',
                f'ensemble npt_mttk temp {temperature} {temperature} iso 0 0',
                #f'mc canonical 100 100 {temperature} {temperature}',
                'dump_exyz 10000 0 1', 
                'run 200000']
        Morph(atoms).gpumd(dirname=dirname, run_in=run_in)
        Temp_Struts += read_xyz(dirname + '/dump.xyz')

MultiMol(Temp_Struts).dump('npt.xyz')
