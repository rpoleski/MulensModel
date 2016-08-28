from astropy import units as u

import MulensModel


# first way
lens = MulensModel.Lens(n_components=3)
lens.s_2 = 1.3
lens.q_2 = 0.001
lens.s_3 = 0.8
lens.q_3 = 0.4
lens.beta_23 = 30. * u.deg

# second way
lens_2 = MulensModel.Lens(n_components=3)
lens_2.set_component(component=2, s=0.5, q=0.4, beta=0.*u.deg)
lens_2.set_component(component=3, s=0.6, q=34., beta=50.*u.deg)

# third way
lens_3 = MulensModel.Lens(n_components=3, s=[1.1, 0.5], q=[0.01, 0.2], beta=[0.0, 30.]*u.deg)

# and some usage:
print(lens.s, lens.q) #lens.s = [lens.s2,lens.s3]
print(lens.s_23)
print(lens.component[2])  # i.e. s2, q2, beta2



