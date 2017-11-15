import MulensModel
"""
Show how a dictionary implementation of ModelParameters would work.
"""

MulensModel.modelparameters.which_parameters('PSPL')
# Returns: 't_0', 'u_0', 't_E'

MulensModel.modelparameters.which_parameters('FSBL')
# Returns: 't_0', 'u_0', 't_E', 'rho', 's', 'q', 'alpha'

MulensModel.modelparameters.which_parameters('BinaryLens')
# Returns: 's', 'q', 'alpha'

MulensModel.modelparameters.which_parameters('BinaryLensOrbitalMotion')
# Returns: 
#    's', 'q', 'alpha', 'dsdt', 'dalphadt' (OPTIONAL: 'z', 'dzdt')

MulensModel.modelparameters.which_parameters()
# Returns: (a more verbose listing)
#    Some common model types:
#        PSPL: 't_0', 'u_0', 't_E'
#        FSPL: 't_0', 'u_0', 't_E', 'rho'
#        PSPL w/ parallax: 't_0', 'u_0', 't_E', 'pi_E_N', 'pi_E_E' (OPTIONAL: 't_0_par')
#        FSBL: 't_0', 'u_0', 't_E', 'rho', 's', 'q', 'alpha'
#        BSPL: 't_0_1', 'u_0_1', 't_0_2', 'u_0_2', 't_E'
#    By Effect:
#        parallax: 'pi_E_N', 'pi_E_E' (OPTIONAL: 't_0_par')
#        xallarap: 'xi_E_N', 'xi_E_E', 'period'
#        finite source: 'rho' or 'rho_1', 'rho_2' (if 2 sources)
#        lens orbital motion: 'dsdt', 'dalphadt' (OPTIONAL: 'z', 'dzdt')
#    Alternative parameterizations:
#        FSPL: 't_0', 't_eff', 't_E', 't_star'
#        FSBL: 't_1', 'u_0', 't_2', 'rho', 's', 'q', 'alpha' (Cassan)

PSPL_params = MulensModel.ModelParameters(
    {'t_0': 2458060., 'u_0': 0.2, 't_E': 30.5})

my_PSPL_model = MulensModel.Model({'t_0': 2458060., 'u_0': 0.2, 't_E': 30.5})
print(my_PSPL_model.parameters)
print(my_PSPL_model.parameters.t_eff)
print(my_PSPL_model.parameters.rho) 
# Returns: AttributeError('rho is not defined for this model.')

import copy

FSPL_params = copy.deepcopy(my_PSPL_model.parameters)
FSPL_params['rho'] = 0.001

my_FSPL_model = MulensModel.Model(FSPL_params)

my_PSPL_model = MulensModel.Model(
    parameters={'t_0': 2458060., 'u_0': 0.2, 't_E': 30.5, 't_0_1': 2458062.})
# Returns: AttributeError('Not a valid combination of parameters. See
#    which_parameters()')
