from astropy.coordinates import SkyCoord
from astropy import units as u


def read_files_from_config(config):
    """
    Read names of photometry files. For each file read also format:
    'mag' or 'flux'.
    """
    section = 'photometry files'
    ids = [t for t in config[section]]
    files = [config.get(section, id_).split() for id_ in ids]
    return files

def read_model_settings(config):
    """
    Read basic parameters of the microlensing models: methods used for
    calculating magnification and coordinates.
    """
    settings = {}
    section = 'model'
    if config.has_option(section, 'methods'):
        methods = config.get(section, 'methods').split()
        for i in range(0, len(methods), 2):
            methods[i] = float(methods[i])
        settings['methods'] = methods
    if config.has_option(section, 'default_method'):
        settings['default_method'] = config.get(section, 'default_method')
    try:
        ra = config.getfloat(section, 'RA')
        dec = config.getfloat(section, 'Dec')
        settings['coords'] = SkyCoord(ra, dec, unit=u.deg)
    except:
        settings['coords'] = None
    return settings

def read_parameters_start(config):
    """
    Read info on starting parameter values. The distributions can be
    'gauss', 'uniform', or 'log-uniform'.
    """
    section = 'EMCEE starting'
    parameters_to_fit = [var for var in config[section]]
    starting = {}
    for param in parameters_to_fit:
        words = config.get(section, param).split()
        starting[param] = [words[0]] + [float(word) for word in words[1:]]
    return (parameters_to_fit, starting)

def read_fix_parameters(config):
    """
    Read parameters that will be kept fixed during fitting process.
    """
    section = 'fixed parameters'
    fixed = {}
    if section not in config:
        return fixed
    for var in config[section]:
        fixed[var] = config.getfloat(section, var)
    return fixed

def read_min_max(config):
    """
    Read minimum and maximum values of parameters used for prior.
    """
    min_values = {}
    section = 'EMCEE min values'
    if section in config.sections():
        for var in config[section]:
            min_values[var] = config.getfloat(section, var)

    max_values = {}
    section = 'EMCEE max values'
    if section in config.sections():
        for var in config[section]:
            max_values[var] = config.getfloat(section, var)
    return (min_values, max_values)

def read_emcee_settings(config):
    """
    Read EMCEE parameters, i.e., number of walkers, steps, and burn-in steps.
    """
    emcee_settings = {}
    section = 'EMCEE'
    for name in ['n_walkers', 'n_steps', 'n_burn']:
        emcee_settings[name] = config.getint(section, name)
    return emcee_settings

def read_other(config):
    """
    Read other options.
    """
    other = {}
    section = 'other'
    if 'print models' in config[section]:
        other['print_models'] = config.getboolean(section, 'print models')
    plot = config.get(section, 'plot time').split()
    other['plot_time'] = [float(t)-2450000. for t in plot]
    return other
