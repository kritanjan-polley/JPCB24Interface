import ForceField as build

# ff = build.ForceField(['swm4-ndp'], polarizable=True, symmetrized=True)
# ff.print('pair-drude_tip4p.lmp')

ff = build.ForceField(['swm4-ndp','O3','polNa','polI','polCl'], polarizable=True, symmetrized=True)
ff.print('pair-drude_tip4p_O3_NaICl.lmp')
