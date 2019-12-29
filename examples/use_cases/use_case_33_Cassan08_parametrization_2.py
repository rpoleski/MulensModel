"""
Use Cassan 2008 parameters for binary lens model.
Translate these parameters to standard ones for a few sets of parameters.
"""
import MulensModel as mm


def nice_print(model):
    """prints selected parameters"""
    print(' t_0 = {:.5f}'.format(model.parameters.t_0))
    print(' u_0 = {:.5f}'.format(model.parameters.u_0))
    print(' t_E = {:.3f}'.format(model.parameters.t_E))
    print(' alpha = {:.2f}\n'.format(model.parameters.alpha))


parameters = {
    's': 1.1, 'q': 0.05,
    'x_caustic_in': 0.4, 'x_caustic_out': 0.65,
    't_caustic_in': 10., 't_caustic_out': 20.}

model = mm.Model(parameters)

print("Original model:")
print(model)
nice_print(model)

model.parameters.x_caustic_in = 0.3
model.parameters.t_caustic_out = 30.
print("Changes:")
print("x_caustic_in -> 0.3")
print("t_caustic_out -> 30.")
print("result in:")
nice_print(model)

model.parameters.s = 1.2
print("Change:")
print("s -> 1.2")
print("results in ...(these calculations can take some time)...")
nice_print(model)

model.parameters.x_caustic_in = 0.25
print("Change:")
print("x_caustic_in -> 0.25")
print("results in:")
nice_print(model)
