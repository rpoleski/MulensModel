import glob
import os
import numpy as np
import sys
import warnings

def get_yaml_name(event, model):
    nazwa_yaml = str(event) + '_' + str(model) + '.yaml'
    return nazwa_yaml

def get_ids(list_of_models):
    ids = []
    for model in list_of_models:
        ids.append(model[-12:-4])
    return ids

def get_parameters_names(ID):
    model_id_short = ID[0:2]
    if model_id_short == 'BO' or model_id_short == 'LK' or model_id_short == 'TS' or model_id_short == 'TX' or model_id_short == 'TO':
        warnings.warn('The transformation of this model category is not accurate!')
        print("This model category: {:}, is not yet included in the transformation.".format(model_id_short))
    model_and_parnames = {'PS': ['u_0', 't_E', 't_0', 'rho'],
                    'PX': ['u_0', 't_E', 't_0', 'rho', 'pi_N', 'pi_E'],
                    'BS': ['t_E', 'FR', 'u_0_1', 'u_0_2', 't_0_1', 't_0_2', 'rho_1'],
                    'BO': ['t_E', 'FR', 'u_0_1', 'u_0_2', 't_0_1', 't_0_2', 'rho_1', 'pi_E_N', 'pi_E_E', 'gamma1', 'gamma2', 'gammaz'],
                    'LS': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0'],
                    'LX': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0', 'pi_E_N', 'pi_E_E'],
                    'LO': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0', 'pi_E_N', 'pi_E_E', 'gamma1', 'gamma2', 'gammaz'],
                    'LK': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0', 'pi_E_N', 'pi_E_E', 'gamma1', 'gamma2', 'gammaz', 'sz_s', 'a_s3d'],
                    'TS': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0', 's2', 'q2', 'beta'],
                    'TX': ['s', 'q', 'u_0', 'alpha', 'rho', 't_E', 't_0', 's2', 'q2', 'beta', 'pi_E_N', 'pi_E_E']}
    for x,y in model_and_parnames.items():
        if x == model_id_short:
            parnames = y
    return parnames

def get_MMformat_for_orbital_motion(params):
    ds_dt = params["s"] * params["gamma1"] * 365.2422
    ds_z_dt = params["gammaz"] * params["s"] * 365.2422
    dalpha_dt = params["gamma2"] * 365.25 *180/np.pi
    params.update({"ds_dt": ds_dt})
    params.update({"ds_z_dt": ds_z_dt})
    params.update({"dalpha_dt": -dalpha_dt})
    new_dict = params.copy()
    new_dict.pop("gamma1")
    new_dict.pop("gamma2")
    new_dict.pop("gammaz")
    return new_dict

def check_t0_format(t0):
    if t0 < 2450000:
        t0 = t0 + 2450000
    else:
        t0 = t0
    return t0    

def get_to_MM_format(params, model):
    if model[0] == "B":
        t01 = params["t_0_1"]
        t01 = check_t0_format(t01)
        params.update({"t_0_1": t01})
        t02 = params["t_0_2"]
        t02 = check_t0_format(t02)
        params.update({"t_0_2": t02})
    else:    
        t0 = params["t_0"]
        t0 = check_t0_format(t0)
        params.update({"t_0": t0})
    if model[0] == "L":
        alpha = params["alpha"] * 180/np.pi
        params.update({"alpha": alpha})
    if model[0:2] == "LO" or model[0:2] == "LK":
        params = get_MMformat_for_orbital_motion(params)
    if model[1] == "X" or model[1] == "O" or model[1] == "K":
        t_0_par = params["t_0"]
        params.update({"t_0_par": t_0_par})
    return params

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Two parameters required:")
        print("1 - direcotry with RTModel results")
        print("2 - 'all' or 'final' - indicates which models should be selected")
        sys.exit(1)

    event = sys.argv[1]
    M = sys.argv[2]
    if M == 'final':
        models = glob.glob(event +'/FinalModels/*')
    elif M == 'all':
        models = glob.glob(event +'/Models/*')
    else:
        raise TypeError("all or final allowed as a second arg")

    models_ids = get_ids(models)
    for i,model in enumerate(models):
        parnames = get_parameters_names(models_ids[i])
        par_dict = dict.fromkeys(parnames, None)
        with open(model) as f:
            lines = f.readlines()
            line = lines[0]
            chunks = line.split(' ')
            values = [float(v) for v in chunks]
            for x, y  in zip(par_dict.keys(), values):
                par_dict.update({x: y})
            par_dict = get_to_MM_format(par_dict, models_ids[i])
            with open(get_yaml_name(event,models_ids[i]), 'x') as file_out:
                print(get_yaml_name(event,models_ids[i]))
                file_out.writelines("Best model:\n")
                file_out.writelines("  Parameters:\n")
                for x,y in par_dict.items():
                    file_out.writelines("    {:}: {:}\n".format(x,y))
