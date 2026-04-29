import glob
import os
import numpy as np
import sys

def get_yaml_name(event, model):
    nazwa_yaml = str(event)[:-1] + '_' + str(model) + '.yaml'
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
        ids.append(model[-12:-6])
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
#def check_values(params):
    #if params["t_0"] < 10:
        #raise TypeError("Something wrong with values")
    #if params["t_E"] < 0:
        #raise TypeError("Something wrong with values")
    #if params["q"] < 0:
        #raise TypeError("Something wrong with values")

if __name__ == '__main__':
    event = sys.argv[1]
    models = glob.glob(event +'/Models/*')

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
            if models_ids[i][0] == "L":    
                alpha = par_dict["alpha"] * 180/np.pi    
                par_dict.update({"alpha": alpha})
            if models_ids[i][0:2] == "LO" or models_ids[i][0:2] == "LK":
                par_dict = get_MMformat_for_orbital_motion(par_dict)
            #print(par_dict.values())
            #check_values(par_dict)
            with open(get_yaml_name(event,models_ids[i]), 'w') as file_out:
                print(get_yaml_name(event,models_ids[i]))
                file_out.writelines("Best model:\n")
                file_out.writelines("  Parameters:\n")
                for x,y in par_dict.items():
                    #print("%s: %f \n" %(x, y))
                    file_out.writelines("    %s: %f \n" %(x, y))
