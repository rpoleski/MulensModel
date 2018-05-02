import MulensModel as mm
import os

raise NotImplementedError('This use case has not been implemented')

data_path = os.path.join(mm.MODULE_PATH, 'data')

ob03235_ogle_data = mm.MulensData(
    data_file=os.path.join(data_path, 'OB03235', 'OB03235_OGLE.tbl.txt'),
    phot_fmt='mag', 
    plot_properties={'marker': 'o', 'size': 5, 'color': 'black', zorder=10})
ob03235_moa_data = mm.MulensData(
    data_file=os.path.join(data_path, 'OB03235', 'OB03235_MOA.tbl.txt'),
    phot_fmt='flux', 
    plot_properties={'marker': 's', 'size': 2, 'color': 'red', zorder=2,
                     show_errorbars=False, show_bad=False})
