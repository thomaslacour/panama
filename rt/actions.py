import csv
import numpy as np

def read(filename, v=0):
    """ Routine to read different kind of file """
    # CONFIG FILE ==================
    if filename=='conf.txt':
        from .config_file import read_config_file
        return read_config_file(filename, v=0)
    else:
        return None
    pass


# def extract_layers_properties(data, lm):
#     """
#     """
#     d = []
#     CL = []
#     CT = []
#     rhoL = []
#     rhoT = []
#     for layer in data['layers'].values():
#         d.append(layer[1])
#         material_name = layer[0]
#         if lm[material_name]['type'] == 'meta':
#             CL.append('meta')
#             CT.append('meta')
#             rhoL.append('meta')
#             rhoT.append('meta')
#         else:
#             CL.append(lm[material_name]['CL'])
#             CT.append(lm[material_name]['CT'])
#             rhoL.append(lm[material_name]['rho'])
#             rhoT.append(lm[material_name]['rho'])
#     if rhoT==[]:
#         return (d, CL, CT, rhoL, rhoL)
#     else:
#         return (d, CL, CT, rhoL, rhoT)




def process_config(data):
    """ Raw data parser
    This function allows communication between computation and raw data.
    """
    from rt import display
    from rt import validation

    try:
        data['layers']
    except KeyError:
        raise ValueError('No layers or wrong syntax in config file !')

    do = data["todo"]
    try:
        do = int(do)
    except:
        pass

    if do == -1:   #validation tests
        validation.test_1layer_NormalIncidence(data)
        validation.test_1layer_Angle(data)
        # fRT = [ abs(u) for u in validation.test_1layer_NormalIncidence(data) ]
        # aRT = [ abs(u) for u in validation.test_1layer_Angle(data) ]
        # csv_field = [*fRT, *aRT]
        # data['todo'] = -1

    elif do == 0:   #RT versus omega
        freq, R, T = display.rt_versus_omega(data)
        return {"freq_MHz":freq, "abs_R":R, "abs_T":T}

    elif do == 1:   #RT versus angle
        theta, R, T = display.rt_versus_theta(data)
        return {"angle_deg":theta, "abs_R":np.abs(R), "abs_T":np.abs(T)}

    elif do == 2:   #RT versus angle and frequency
        pass

    elif do == 3:   #transfert matrix versus frequency
        display.transfert_matrix_versus_omega(data)

    elif do == 4:   #scattering matrix versus frequency
        display.scattering_matrix_versus_omega(data)

    elif type(do)!=int and do[0:4] == 'PROP': #effective properties of 1 layer
        layer_number = do[4:]
        layer_prop = data['layers'][layer_number]
        from material import display as matplt
        f = np.linspace(data['f_min'], data['f_max'], data['f_num'])
        (cLeff,alphaLeff,cTeff,alphaTeff) = matplt.prop(layer_prop, f)
        csv_field = {"freq":f, "cL":cLeff, "alphaL":alphaLeff, "cT":cTeff, "alphaT":alphaTeff}

    else:
        print("Unable to find the action to execute.")

    return




def save_as_csv(things, csvfile):
    """
    Pandas module lets you save a "things" to a csv file.
    If the things is a dict, the row header are dict keys.
    Otherwise, list is save without header.
    /!\ remark: datas are formated with 5 digits.
    """
    try:
        import pandas as pd
    except ModuleNotFoundError:
        return

    if type(things)==dict:
        list_of_thing = list(things.values())
        header = [*things]
    elif type(things)==list:
        list_of_thing = things
        header=['']*len(things)
    elif things==None:
        return
    my_df = pd.DataFrame(list_of_thing)
    my_df = my_df.T
    my_df.to_csv(csvfile, header=header, index=False, float_format='%.5e' )
