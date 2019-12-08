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

# There are 3 approaches:
# First approach - a separate function for each degeneracy:
model_2 = model_1.predict_degenerate_jerk_model()
model_3 = model_2.predict_degenerate_u_0_model()

# Second approach - have a single function that deals with all degeneracies:
model_2 = model_1.predict_degenerate_model(jerk=True)
# and this allows combining multiple degeneracies:
model_3 = model_1.predict_degenerate_model(jerk=True, u_0_sign=True)
# other parameters we could have here are e.g.,
# "close_wide" (would work only for binary lens events),
# "ecliptic"...

# Third approach - functions predicting degeneracies return a dict with
# parameters to be changed and if the user wants to use them in
# a Model instance, then they should copy the Model instance (which is not
# implemented at this point) and change specific parameters.
parameters_2 = model_1.predict_degenerate_jerk_model()
parameters_3 = {'u_0': -model_1.parameters.u_0}
model_2 = model_1.copy()
# This copies not only parameters but coordinates etc.
model_3 = model_1.copy()
for (key, value) in parameters_2:
    setattr(model_2.parameters, key, value)
for (key, value) in parameters_3:
    setattr(model_3.parameters, key, value)

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
