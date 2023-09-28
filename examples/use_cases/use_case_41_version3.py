"""
Use case for Version 3 model sub-classes.

For convenience, it would be a good idea to have both full names and a
short-hand version of the model classes.
"""
import MulensModel as mm

par_params = {'pi_E_N': 0.1, 'pi_E_E': 0.3}
coords = '18:00:00 -30:00:00'

# point lens models:
pspl_params = {'t_0': 0., 'u_0': 0.001, 't_E': 1}
pspl_params_w_par = {**pspl_params, **par_params}
fspl_params = {**pspl_params, 'rho': 0.002}
fspl_params_w_par = {**fspl_params, **par_params}

# binary lens models:
planet_params = {'s': 0.95, 'q': 0.001, 'alpha': 90.}
planet_ps_params = {**pspl_params, **planet_params}
planet_ps_params_w_par = {**planet_ps_params, **par_params}
planet_params = {**planet_ps_params, 'rho': 0.0001}
planet_params_w_par = {**planet_params, **par_params}

# binary source models:
# point lens
bspl_params = {'t_0_1': 0., 'u_0_1': 0.001, 't_0_2': -0.3, 'u_0_2': 0.1,
               't_E': 1}
bspl_fs_params_1 = {**bspl_params, 'rho_1': 0.002}
bspl_fs_params_2 = {**bspl_params, 'rho_1': 0.002, 'rho_2': 0.3}
# planet
bsbl_params = {**bspl_params, **planet_params}
bsbl_fs_params_1 = {**bspl_fs_params_1, **planet_params}
bsbl_fs_params_2 = {**bspl_fs_params_2, **planet_params}

list_of_params = [pspl_params, fspl_params, planet_ps_params, planet_params,
                  pspl_params_w_par, fspl_params_w_par, planet_ps_params_w_par,
                  planet_params_w_par,
                  bspl_params, bspl_fs_params_1, bspl_fs_params_2,
                  bsbl_params, bsbl_fs_params_1, bsbl_fs_params_2]

# Basic Model Definition, same as pre
for params in list_of_params:
    if 'pi_E_N' in params.keys():
        model = mm.Model(params, coords=coords)
    else:
        model = mm.Model(params)
    
# Specific models: Short name
pspl_model_1 = mm.PSPLModel(pspl_params)
pspl_model_2 = mm.PSPLModel(pspl_params_w_par, coords=coords)
# Alternative:
pspl_model_1b = mm.PS1LModel(pspl_params)

fspl_model_1 = mm.FSPLModel(fspl_params)
fspl_model_2 = mm.FSPLModel(fspl_params_w_par, coords=coords)
# Alternative:
fspl_model_1b = mm.FS1LModel(pspl_params)

bspl_model_1 = mm.BSPLModel(bspl_params)
bspl_fs_model_2 = mm.BSPLModel(bspl_fs_params_1, coords=coords)
bspl_fs_model_3 = mm.BSPLModel(bspl_fs_params_2, coords=coords)
# Is this the right behavior for binary lenses?
# Alternatives:
bspl_model_1b = mm.BS1LModel(pspl_params)
bspl_model_2b = mm.BSFS1LModel(pspl_params)  # see Long names below.
bspl_model_2c = mm.BFS1LModel(pspl_params)

planet_ps_model_1 = mm.PSBLModel(planet_ps_params)
planet_ps_model_2 = mm.PSBLModel(planet_ps_params_w_par, coords=coords)
# Alternative:
planet_ps_model_1b = mm.PS2LModel(planet_ps_params)

planet_model_1 = mm.FSBLModel(planet_params)
planet_model_2 = mm.FSBLModel(planet_params_w_par, coords=coords)
# Alternative:
planet_model_1b = mm.FS2LModel(planet_ps_params)

bsbl_model_1 = mm.BSBLModel(bsbl_params)
bsbl_fs_model_2 = mm.BSBLModel(bsbl_fs_params_1, coords=coords)
bsbl_fs_model_3 = mm.BSBLModel(bsbl_fs_params_2, coords=coords)
# Is this the right behavior for binary lenses?

# Specific models: Long name
pspl_model_1 = mm.PointSourcePointLensModel(pspl_params)
pspl_model_2 = mm.PointSourcePointLensModel(pspl_params_w_par, coords=coords)

fspl_model_1 = mm.FiniteSourcePointLensModel(fspl_params)
fspl_model_2 = mm.FiniteSourcePointLensModel(fspl_params_w_par, coords=coords)

bspl_model_1 = mm.BinarySourcePointLensModel(bspl_params)
bspl_fs_model_2 = mm.BinarySourceFiniteSourcePointLensModel(
    bspl_fs_params_1, coords=coords)
bspl_fs_model_3 = mm.BinarySourceDualFiniteSourcePointLensModel(
    bspl_fs_params_2, coords=coords)
# Is this the right behavior for binary lenses? This version has 3 different
# model types for Binary Sources rather than just one. They would be
# abbreviated: BSPL, BSFSPL, BS2FSPL.

planet_ps_model_1 = mm.PointSourceBinaryLensModel(planet_ps_params)
planet_ps_model_2 = mm.PointSourceBinaryLensModel(planet_ps_params_w_par, coords=coords)

planet_model_1 = mm.FiniteSourceBinaryLensModel(planet_params)
planet_model_2 = mm.FiniteSourceBinaryLensModel(planet_params_w_par, coords=coords)

bsbl_model_1 = mm.BinarySourceBinaryLensModel(bspl_params)
bsbl_fs_model_2 = mm.BinarySourceFiniteSourceBinaryLensModel(
    bspl_fs_params_1, coords=coords)
bsbl_fs_model_3 = mm.BinarySourceDualFiniteSourceBinaryLensModel(
    bspl_fs_params_2, coords=coords)
