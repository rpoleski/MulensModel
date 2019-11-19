"""
Use cases showing how MulensModel can predict parameters of degenerate models.

Here I use the event ob140939 as an example. I know that jerk parallax,
in general, works for small u_0 events, but I don't have a better
dataset at hand.
"""
import os

import MulensModel as mm


parameters = {
    't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
    'pi_E_N': -0.248, 'pi_E_E': 0.234}
coords = "18:04:45.71 -26:59:15.2"

model_1 = mm.Model(parameters, coords=coords)

# There are 2 approaches: first - a separate function for each degeneracy:
model_2 = model_1.find_degenerate_jerk_model()
model_3 = model_2.find_degenerate_u_0_model()
# second - have a single function that deals with all degeneracies:
model_2 = model_1.find_degenerate_model(jerk=True)
# and this allows combining multiple degeneracies:
model_3 = model_1.find_degenerate_model(jerk=True, u_0_sign=True)
# other parameters we could have here are e.g.,
# "close_wide" (would work only for binary lens events),
# "ecliptic"...

print("the 3 models are:")
print(model_1)
print(model_2)
print(model_3)

# Get the chi2 values for these models:
file_name = os.join(mm.DATA_PATH, "photometry_files", "OB140939",
                    "ob140939_OGLE.dat")
data = mm.MulensData(file_name=file_name)
models = [model_1, model_2, model_3]
for model in models:
    event = mm.Event(datasets=data, model=model)
    print(event.get_chi2())
