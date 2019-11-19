"""
Calculate the Einstein ring for a lens at a variety of distances with
a variety of masses and output the results to a file.

Output:
    thetaE.tab - Angular Einstein ring radius in milliarcseconds

    rE.tab - Physical Einstein ring radius in AU (Lens plane)

    rEtilde.tab - Projection of Einstein ring radius in AU (Observer plane)

    tE.tab - Einstein timescale in days for a given mu_rel

"""
import numpy as np
import astropy.units as u

import MulensModel as mm

# Define the source and lens parameters
lens_dist = np.arange(1., 8., 1.)
lens_mass = [10., 1., 0.3, 0.1, 0.01, 0.001]
source = mm.Source(distance=8.)
mu_rel = 4.

names = ['thetaE.tab', 'rE.tab', 'tE.tab', 'rEtilde.tab']

# Open the output files and print headers
file_thetaE = open(names[0], 'w')
file_rE = open(names[1], 'w')
file_tE = open(names[2], 'w')
file_rEtilde = open(names[3], 'w')
file_thetaE.write('# Angular Einstein ring radius in milliarcseconds.\n')
file_rE.write('# Physical Einstein ring radius in AU (Lens plane).\n')
file_rEtilde.write(
    '# Projection of Einstein ring radius in AU (Observer plane).\n')
file_tE.write(
    '# Einstein timescale in days for mu_rel = {0:4.1f} mas/yr.\n'.format(
        mu_rel))

file_list = [file_thetaE, file_rE, file_rEtilde, file_tE]

# Print column headings
for file_ in file_list:
    file_.write('# columns = distance in kpc\n# rows = mass in solMass\n')
    file_.write('{0:6}'.format(' '))
    for dist in lens_dist:
        file_.write('{0:6.2f}'.format(dist))
    file_.write('\n')

# Calculate the Einstein Ring for a given lens mass...
for mass in lens_mass:
    for file_ in file_list:
        file_.write('{0:6.4}'.format(mass))

    # ...and a given distance
    for dist in lens_dist:
        # Setup the microlensing system
        system = mm.MulensSystem(
            lens=mm.Lens(mass=mass, distance=dist), source=source,
            mu_rel=mu_rel)

        # Output the calculated Einstein radius to a file
        file_thetaE.write('{0:6.2f}'.format(system.theta_E.value))
        file_rE.write('{0:6.2f}'.format(system.r_E.value))
        file_rEtilde.write('{0:6.2f}'.format(system.r_E_tilde.value))
        file_tE.write('{0:6.1f}'.format(system.t_E.value))

    # Add new line characters to files
    for file_ in file_list:
        file_.write('\n')

# Close the files
for file_ in file_list:
    file_.close()

print("Congratulations!")
msg = "Have a look at following files: {:}, {:}, {:}, and {:}"
print(msg.format(*names))
