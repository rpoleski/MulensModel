"""
Extract Roman positions for Data Challenge 2026 in a format known to MulensModel.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d


def get_Earth_xyz(times):
    """
    For given time vector extract positions of L2 relative to Solar System Barycenter.
    """
    padding = 0.2 # in days
    step = "72m" # 0.05 day
    keys = ['x', 'y', 'z']

    begin = int(min(times)/padding) * padding
    end = int(max(times)/padding+1) * padding
    t0 = Time(begin, format='jd')
    t1 = Time(end, format='jd')
    epochs = {'start': t0.iso, 'stop': t1.iso, 'step': step}
    # 399 - Earth
    # 500@0 - Solar System Barycenter
    horizon = Horizons(id='399', location='500@0', epochs=epochs)
    vectors = horizon.vectors()
    kwargs = {'x': vectors['datetime_jd'], 'kind': 'cubic'}
    interp = {key: interp1d(y=vectors[key], **kwargs) for key in keys}
    out = {key: interp[key](times) for key in keys}
    return out


def get_masked_bjd_xyz(dataframe, bjd_begin, bjd_end):
    """
    Mask a pandas data frame based on bjd, select unique values, and extract corresponding x,y,z.
    """
    keys = ['x', 'y', 'z']
    root = "obs_"

    bjd = dataframe['bjd']
    mask = (bjd > bjd_begin) & (bjd <= bjd_end)
    bjd_np = bjd[mask].to_numpy()
    indexes = np.unique(bjd_np, return_index=True)[1]
    bjd_masked = bjd_np[indexes]
    xyz = {key: df[root+key][mask].to_numpy()[indexes] for key in keys}
    return (bjd_masked, xyz)


def transform_GeoEcl_to_GeoEqu(relative):
    """
    Go from:
    geocentric ecliptic cartesian
    to:
    geocentric equatorial cartesian

    Input is a dict with 'x', 'y', 'z' keywords and AU units
    """
    skycoord = SkyCoord(**relative, frame='geocentrictrueecliptic', representation_type='cartesian', unit='au')
    skycoord = skycoord.gcrs
    skycoord.representation_type='cartesian'
    return (skycoord.x.value, skycoord.y.value, skycoord.z.value)


def save_to_file(file_name, bjd, out):
    """save 4-column text file"""
    fmt = "{:.5f} {:.6f} {:.6f} {:.6f}\n"
    with open(file_name, 'x') as out_file:
        for bjd_xyz in zip(bjd, *out):
            out_file.write(fmt.format(*bjd_xyz))


if __name__ == '__main__':
    data_file = "RMDC26_Beginner_Tier_test.parquet"
    out_files_format = "RMDC26_ephemeris_season_{:}.dat"

    borders = [-1000, 1600, 1800, 2000, 2850, 3000, 9000]  # These values allows breaking the data into observing seasons.
    borders = 2460000 + np.array(borders)

    df = pd.read_parquet(data_file)
    df.sort_values('bjd', inplace=True)

    # The loop below goes over seasons.
    for (i, (b0, b1)) in enumerate(zip(borders[:-1], borders[1:])):
        (bjd_masked, xyz) = get_masked_bjd_xyz(df, b0, b1)

        L2 = get_Earth_xyz(bjd_masked)

        diff = {key: xyz[key] - L2[key] for key in ['x', 'y', 'z']}

        out = transform_GeoEcl_to_GeoEqu(diff)

        file_out = out_files_format.format(i+1)
        save_to_file(file_out, bjd_masked, out)

