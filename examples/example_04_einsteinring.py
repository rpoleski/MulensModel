import numpy as np
import astropy.units as u
import MulensModel as mm

"""
Calculate the Einstein ring for a lens at a variety of distances with
a variety of masses and output the results to a file.

Output:
 thetaE.tab - Angular Einstein ring radius in milliarcseconds
 rE.tab - Physical Einstein ring radius in AU (Lens plane)
 rEtilde.tab - Projection of Einstein ring radius in AU (Observer plane)
 tE.tab - Einstein timescale in days for a given mu_rel
"""

#Define the source and lens parameters
lens_dist = np.arange(1., 8., 1.) 
lens_mass = [10., 1., 0.3, 0.1, 0.01, 0.001]
source = mm.Source(distance=8.)
mu_rel = 4.

#Open the output files and print headers
file_thetaE = open('thetaE.tab', 'w')
file_rE = open('rE.tab', 'w')
file_tE = open('tE.tab', 'w')
file_rEtilde = open('rEtilde.tab', 'w')
file_thetaE.write('# Angular Einstein ring radius in milliarcseconds.\n')
file_rE.write('# Physical Einstein ring radius in AU (Lens plane).\n')
file_rEtilde.write('# Projection of Einstein ring radius in AU (Observer plane).\n')
file_tE.write('# Einstein timescale in days for mu_rel = {0:4.1f} mas/yr.\n'.format(mu_rel))

file_list = [file_thetaE, file_rE, file_rEtilde, file_tE]

#Print column headings
for file in file_list:
    file.write('# columns = distance in kpc\n# rows = mass in solMass\n')
    file.write('{0:6}'.format(' '))
    for dist in lens_dist:
        file.write('{0:6.2f}'.format(dist))
    file.write('\n')

#Calculate the Einstein Ring for a given lens mass...
for mass in lens_mass:
    for file in file_list:
        file.write('{0:6.4}'.format(mass))

    #...and a given distance
    for dist in lens_dist:
        #Setup the microlensing system
        system = mm.MulensSystem(lens=mm.Lens(mass=mass,distance=dist),
                                 source=source,mu_rel=mu_rel)

        #Output the calculated Einstein radius to a file
        file_thetaE.write('{0:6.2f}'.format(system.theta_E.value))
        file_rE.write('{0:6.2f}'.format(system.r_E.value))
        file_rEtilde.write('{0:6.2f}'.format(system.r_E_tilde.value))
        file_tE.write('{0:6.1f}'.format(system.t_E.value))

    #Add new line characters to files
    for file in file_list:
        file.write('\n')

#Close the files
for file in file_list:
    file.close()
