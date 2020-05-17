import numpy as np
import copy as cp 

eps = 1e-9

def space_parameters_ND_meshed(eTotal, phiEff, eMax, eMin, eStep,
                                               phiMax, phiMin, phiStep,
                                               eMaxLast, eMinLast,
                                               phiMaxLast, phiMinLast):
    """
    SPACE_PARAMETERS_ND: non monotonous parametrization of the porous multilayer.

    Define all the possible combinations of thickness and concentration
    for a given set of constraints (mins, maxs, steps, and totals/effectives).

    "eTotal" is sum of all thickness and "phiEff" is the effective concentration
    of the multilayer. For each layer, thickness and concentration are bounded
    with a "Max" and a "Min" (so eMax, eMin, eStep, phiMax, phiStep are vectors).

    The last layer is treated separetely becaus its "e" and its "phi" depends
    on the other "front" layers.

    The returned "listOfParam" is a tuple with (E0, PHI0, E1, PHI1, ... ), where
    each element is ndarray of the same size/shape.

    The ndarray "param_valid" is a boolean array with valid parameters satisfying
    the constraints. (size of ndarray param)
    """

    #Generating the stepped data for each dimension
    n = len(eMax) #n is the number of front layers (total number of layers = n+1)
    e = (np.arange(eMin[i], eMax[i]+eps, eStep[i]) for i in range(n))
    phi = (np.arange(phiMin[i], phiMax[i]+eps, phiStep[i]) for i in range(n))
    #Generating the Meshgrids
    meshed_data = np.meshgrid(*e, *phi, indexing='ij')
    E_last = eTotal
    for i in range(n):
        E_last = E_last - meshed_data[i]
    #Preventing 0 values from E_last while dividing by E_last (this case will be discarded later anyway)
    E_last[E_last==0] = -1.0
    PHI_last = phiEff*eTotal
    for i in range(n):
        PHI_last = PHI_last - meshed_data[i+n]*meshed_data[i]
    PHI_last = PHI_last/E_last

    #Logical selections for satisfying constraints
    if eMinLast > 0.0:
        eValid = np.logical_and(E_last >= eMinLast, E_last <= min(eTotal,eMaxLast))
    else:
        eValid = np.logical_and(E_last > 0.0, E_last <= min(eTotal,eMaxLast))
    phiValid = np.logical_and(PHI_last >= phiMinLast, PHI_last <= phiMaxLast)
    param_valid = np.logical_and(eValid, phiValid)

    #Deselecting all non-zero values of phi0, phi1, phi2... when e0, e1, e2... = 0
    Full_indices = [slice(None,None) for j in range(2*n)] #all indices range from 0 to their last value
    for i in range(n):
        if eMin[0] == 0.0: #the index 0 represents the least value of ith e
            indices_Layer_i = cp.copy(Full_indices)
            #Selection of the first value (0) of the ith thickness dimension
            indices_Layer_i[i] = slice(0,1)
            #Selection of all but the first values of ith phi
            indices_Layer_i[i+n] = slice(1,None)
            #Deselecting all non-zero values of ith phi when ith e = 0
            indices_Layer_i = tuple(indices_Layer_i)
            param_valid[indices_Layer_i] = False

    listOfParamMesh = tuple(meshed_data[j//2] if j%2 == 0 else meshed_data[n+(j-1)//2] for j in range(2*n)) + (E_last,PHI_last)

    return (listOfParamMesh, param_valid)



def space_parameters_ND(eTotal, phiEff, eMax, eMin, eStep,
                                        phiMax, phiMin, phiStep,
                                        eMaxLast, eMinLast,
                                        phiMaxLast, phiMinLast):
    """
    Conversion of the list of data_meshed into a list of tuple of valid parameters.
    The output format is

        [ (e0_1, phi0_1, e1_1, ...), (e0_2, phi0_2, e1_2, ...), ... ]

    """
    (data_meshed, param_valid) = space_parameters_ND_meshed(eTotal, phiEff,
            eMax, eMin, eStep,
            phiMax, phiMin, phiStep,
            eMaxLast, eMinLast,
            phiMaxLast, phiMinLast)
    list_of_valid = []
    data_valid = tuple(item[param_valid] for item in data_meshed)
    for item in zip(*data_valid):
        list_of_valid.append(item)
    return list_of_valid



def space_parameters_mono(eTotal, phiEff, phiMax, phiMin,
                                          qmax, qmin, qstep
                                          ):
    """
    SPACE_PARAMETERS_MONO: monotonous parametrization of a porous multilayer.

        phi(x) = A*x^q     --->    phi[n] = phiEff/N^q * (n^(q+1)-(n-1)^(q+1))

    The thickness of each layer is forced to be constant. Only the concentration in
    the layers are allowed to change and:  max(phiMin,0) < phi < phiMax
    The output format is:

        [ (e, phi0_q0, e, phi1_q0, ...), (e, phi0_q1, e, phi1_q1, ...), ... ]

    """
    list_of_param = []
    N = len(phiMin)
    e = eTotal/N
    q_list = np.arange(qmin, qmax, qstep)
    n = np.arange(N) + 1
    for q in q_list:
        phin = phiEff/N**q * ( n**(q+1) - (n-1)**(q+1) )
        param = tuple( phin[j//2] if j%2!=0 else e for j in range(2*N))
        list_of_param.append(param)
    return list_of_param



if __name__ == "__main__":

    eTotal = 10.0
    phiEff = 10.0

    nb_front_layers = 2

    eMax = (eTotal,)*nb_front_layers
    eMin = (0.0,)*nb_front_layers
    eStep = (eTotal/5,)*nb_front_layers
    phiMax = (30.0,)*nb_front_layers
    phiMin = (0.0,)*nb_front_layers
    phiStep = (phiMax[0]/4,)*nb_front_layers
    eMaxLast, eMinLast = eTotal, 0.0
    phiMaxLast, phiMinLast = 10.0, 0.0

    list_of_param = space_parameters_ND(eTotal, phiEff, eMax, eMin, eStep,
                                                        phiMax, phiMin, phiStep,
                                                        eMaxLast, eMinLast,
                                                        phiMaxLast, phiMinLast)

    # for u in list_of_param: print(u)

    list_of_param = space_parameters_mono(eTotal, phiEff, (30,)*10, (0,)*10, 2, 1, 0.1)
