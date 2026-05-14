import glob
import os
import numpy as np
import sys

def get_yaml_name(event, model):
    nazwa_yaml = str(event) + '_' + str(model) + '.yaml'
    return nazwa_yaml

def get_index(list_of_models):
    modelcodes= ['PS','PX','BS','BO','LS','LX','LO','LK','TS','TX']
    idx = []
    for model in list_of_models:
        model_id_short = model[-12:-10] #from modelcodes
        idx.append(modelcodes.index(model_id_short)) #index to use for parnames
    return idx

def get_ids(list_of_models):
    ids = []
    for model in list_of_models:
        ids.append(model[-12:-4])
    return ids
def get_parameters_names(index):
    parnames = [['u_0','t_E','t_0','rho'],
                    ['u_0','t_E','t_0','rho','pi_N','pi_E'],
                    ['t_E','FR','u_0_1','u_0_2','t_0_1','t_0_2','rho_1'],
                    ['t_E','FR','u_0_1','u_0_2','t_0_1','t_0_2','rho_1','pi_E_N','pi_E_E','gamma1','gamma2','gammaz'],
                    ['s','q','u_0','alpha','rho','t_E','t_0'],
                    ['s','q','u_0','alpha','rho','t_E','t_0','pi_E_N','pi_E_E'],
                    ['s','q','u_0','alpha','rho','t_E','t_0','pi_E_N','pi_E_E','gamma1','gamma2','gammaz'],
                    ['s','q','u_0','alpha','rho','t_E','t_0','pi_E_N','pi_E_E','gamma1','gamma2','gammaz','sz_s','a_s3d'],
                    ['s','q','u_0','alpha','rho','t_E','t_0','s2','q2','beta'],
                    ['s','q','u_0','alpha','rho','t_E','t_0','s2','q2','beta','pi_E_N','pi_E_E']]
    return parnames[index]

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

def get_to_MM_format(params, model):
    t0 = params["t_0"]
    if t0 < 2450000:
        t0 = t0 + 2450000
    params.update({"t_0": t0})
    if model[0] == "L":
        alpha = params["alpha"] * 180/np.pi
        params.update({"alpha": alpha})
    if model[0:2] == "LO" or model[0:2] == "LK":
        params = get_MMformat_for_orbital_motion(params)
    if model[1] == "X":
        t_0_par = par_dict["t_0"]
        par_dict.update({"t_0_par": t_0_par})
    return params

if __name__ == '__main__':
    event = sys.argv[1]
    M = sys.argv[2]
    if M == 'final':
        models = glob.glob(event +'/FinalModels/*')
    elif M == 'all':
        models = glob.glob(event +'/Models/*')
    else:
        raise TypeError("all or final allowed as a second arg")

    models_idxs = get_index(models)
    models_ids = get_ids(models)
    for i,model in enumerate(models):
        parnames = get_parameters_names(models_idxs[i])
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
