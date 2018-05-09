import MulensModel as mm

import os
import glob
import matplotlib.pyplot as pl

raise NotImplementedError('This use case has not been implemented')

data_path = os.path.join(mm.MODULE_PATH, 'data', 'photometry_files')

# Basic: Two datasets with specified data properties
ob03235_ogle_data = mm.MulensData(
    data_file=os.path.join(data_path, 'OB03235', 'OB03235_OGLE.tbl.txt'),
    phot_fmt='mag', 
    plot_properties={'size': 5, 'color': 'black', 'zorder'=10})
ob03235_moa_data = mm.MulensData(
    data_file=os.path.join(data_path, 'OB03235', 'OB03235_MOA.tbl.txt'),
    phot_fmt='flux', 
    plot_properties={'marker': 's', 'size': 2, 'color': 'red', 'zorder'=2,
                     'show_errorbars'=False, 'show_bad'=False})

# Set plot_properties for many datasets
def set_plot_properties(filename):
    """
    Get "standard" plot properties based on the filename.

    :param filename:
    :return plot_properties:
    """

    plot_properties = {}
    if 'OGLE' in filename:
        plot_properties['color']  = 'black'
        plot_properties['zorder'] = 10
    elif 'MOA' in filename:
        plot_properties['color'] = 'red'
        plot_properties['zorder'] = 2
        plot_properties['show_errorbars'] = False
    else:
        return None

    return plot_properties


file_list = glob.glob(os.path.join(data_path, 'MB08310', '*'))
datasets = []
for file in file_list:
    plot_properties = set_plot_properties(file)
    if plot_properties == None:
        plot_properties = {}

    plot_properties['label'] = file.split('_', maxsplit=2)[0]
    datasets.append(mm.MulensData(file_name=file), plot_properties=plot_properties)

t_0 = 2454656.39975
u_0 = 0.00300
t_E = 11.14
t_star = 0.05487
model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 't_star': t_star})

event = mm.Event(datasets=datasets, model=model)
event.plot_data()
pl.show()