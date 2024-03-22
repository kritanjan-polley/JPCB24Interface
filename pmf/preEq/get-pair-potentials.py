import ForceField as build

# for water + NaI + NaCl
ff = build.ForceField(['swm4-ndp','O3','polNa','polI','polCl'], polarizable=True, symmetrized=True)
ff.print('pair-drude_tip4p_O3_NaICl.lmp')

# for pure water 
# ff = build.ForceField(['swm4-ndp','O3','polNa'], polarizable=True, symmetrized=True)
# ff.print('pair-drude_tip4p_O3_water.lmp')

# for  water + NaI
# ff = build.ForceField(['swm4-ndp','O3','polNa','polI'], polarizable=True, symmetrized=True)
# ff.print('pair-drude_tip4p_O3_NaI.lmp')
