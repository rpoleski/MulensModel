"""
Shows how a dictionary implementation of ModelParameters would work.
"""
import MulensModel as mm


raise NotImplementedError('which_parameters() was removed in version 3.0.0')

mm.modelparameters.which_parameters('PSPL')
# Returns: 't_0', 'u_0', 't_E'

mm.modelparameters.which_parameters('FSBL')
# Returns: 't_0', 'u_0', 't_E', 'rho', 's', 'q', 'alpha'

mm.modelparameters.which_parameters('BinaryLens')
# Returns: 's', 'q', 'alpha'

mm.modelparameters.which_parameters('BinaryLensOrbitalMotion')
# Returns:
#    's', 'q', 'alpha', 'dsdt', 'dalphadt' (OPTIONAL: 'z', 'dzdt')

mm.modelparameters.which_parameters()
# Returns: (a more verbose listing)
#    Some common model types:
#        PSPL: 't_0', 'u_0', 't_E'
#        FSPL: 't_0', 'u_0', 't_E', 'rho'
#        PSPL w/ parallax: 't_0', 'u_0', 't_E', 'pi_E_N', 'pi_E_E'
#            (OPTIONAL: 't_0_par')
#        FSBL: 't_0', 'u_0', 't_E', 'rho', 's', 'q', 'alpha'
#        BSPL: 't_0_1', 'u_0_1', 't_0_2', 'u_0_2', 't_E'
#    By Effect:
#        parallax: 'pi_E_N', 'pi_E_E' OR 'pi_E' (OPTIONAL: 't_0_par')
#        xallarap: 'xi_E_N', 'xi_E_E', 'period'
#        finite source: 'rho' (1 source) OR 'rho_1', 'rho_2' (if 2 sources)
#        lens orbital motion: 'dsdt', 'dalphadt' (OPTIONAL: 'z', 'dzdt')
#    Alternative parameterizations:
#        any two of 'u_0', 't_E', 't_eff' (t_eff = u_0 * t_E)
#        any two of 't_E', 'rho', 't_star' (t_star = rho * t_E)
#        FSBL: 't_1', 'u_0', 't_2', 'rho', 's', 'q', 'alpha' (Cassan)

PSPL_params = mm.ModelParameters(
    {'t_0': 2458060., 'u_0': 0.2, 't_E': 30.5})

my_PSPL_model = mm.Model({'t_0': 2458060., 'u_0': 0.2, 't_E': 30.5})
print(my_PSPL_model.parameters)
print(my_PSPL_model.parameters.t_eff)
print('rho:', my_PSPL_model.parameters.rho)
# Returns: None

FSPL_params = dict(my_PSPL_model.parameters.parameters)
FSPL_params['rho'] = 0.001

my_FSPL_model = mm.Model(FSPL_params)

try:
    my_PSPL_model = mm.Model(
        parameters={
            't_0': 2458060., 'u_0': 0.2, 't_E': 30.5, 't_0_1': 2458062.})
except ValueError as msg:
    print('ValueError: ', msg)
# Returns: ValueError('Not a valid combination of parameters. See
#    which_parameters()')
