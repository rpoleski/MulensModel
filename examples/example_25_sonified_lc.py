"""
Example usage of MulensModel to sonify light curve data file.
data file phot_ob08092_O4.dat.
"""

import MulensModel as mm
from astronify.series import SoniSeries
from astropy.table import Table
import os.path

print('Not tested due to problems installing astronify. -- JCY')

# Read in the data file
SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "photometry_files", "OB08092",
    "phot_ob08092_O4.dat")
data = mm.MulensData(file_name=SAMPLE_FILE_01)

# Sonify
data_table = Table({"time": data.time,
                    "flux": data.flux})
data_soni = SoniSeries(data_table)
data_soni.note_spacing = 0.2
data_soni.sonify()
data_soni.play()
