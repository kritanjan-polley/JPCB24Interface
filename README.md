# JPCB24a

## pmf 
Contains scripts to compute the free energy profiles of pulling an ozone through an air-water interface with umbralla sampling. 
### preEq 
Makes the initial configuration of the slab geometry and then adds drude ocillator to polarizable atoms (at 0K).
### equilibration
Takes the initial polarizable configuration from preEq folder, heats it up slowly and then equilibrates.
### production
Takes equilibrated data files from equilibration folder and runs final production. The final production data files (lammps output format) will have name data.${location}.prod.pol

## shooting 
Requires input data files from the pmf folder (not provided here). The shooting trajectories will generate data files containing 3 columns, time (fs), location (A), and velocity (A/fs)
