"""
Use cases showing how MulensModel can predict parameters of degenerate models.

Here I use the event ob140939 as an example. I know that jerk parallax,
in general, works for small u_0 events, but I don't have a better
dataset at hand.
"""
import os

import MulensModel as mm


raise NotImplementedError('predict_degenerate_jerk_model not implemented')

parameters = {
    't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
    'pi_E_N': -0.248, 'pi_E_E': 0.234}
coords = "18:04:45.71 -26:59:15.2"

model_1 = mm.Model(parameters, coords=coords)

# Proposed solution - a separate function for each degeneracy:
model_2 = model_1.predict_degenerate_jerk_model()
# There are two other degeneracies that can be implemented this way
# (i.e., discrete and requiring non-trivial calculations):
# 1) satellite - Gould 1994 (1994ApJ...421L..75G) and
#    Refsdal 1966 (1966MNRAS.134..315R),
# 2) rho - Chung et al. 2017 (2017ApJ...838..154C).
# References for jerk degeneracy:
# Gould 2004 (2004ApJ...606..319G) and Park et al. 2004 (2004ApJ...609..166P).

print("the 2 models are:")
print(model_1)
print(model_2)

# Get the chi2 values for these models:
file_name = os.join(mm.DATA_PATH, "photometry_files", "OB140939",
                    "ob140939_OGLE.dat")
data = mm.MulensData(file_name=file_name)
for model in [model_1, model_2]:
    event = mm.Event(datasets=data, model=model)
    print(event.get_chi2())
