import numpy as np
import copy as cp 

eps = 1e-9

def space_parameters_ND(eTotal, phiEff, eMax, eMin, eStep,
                                        phiMax, phiMin, phiStep,
                                        eMaxLast, eMinLast,
                                        phiMaxLast, phiMinLast):
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

    listOfParam = tuple(meshed_data[j//2] if j%2 == 0 else meshed_data[n+(j-1)//2] for j in range(2*n)) + (E_last,PHI_last)

    return (listOfParam, param_valid)







if __name__ == "__main__":

    eTotal = 10.0
    phiEff = 5.0

    nb_front_layers = 2

    eMax = (eTotal,)*nb_front_layers
    eMin = (0.0,)*nb_front_layers
    eStep = (2.0,)*nb_front_layers
    phiMax = (10.0,)*nb_front_layers
    phiMin = (0.0,)*nb_front_layers
    phiStep = (2.0,)*nb_front_layers
    eMaxLast, eMinLast = eTotal, 0.0
    phiMaxLast, phiMinLast = 10.0, 0.0
    print('ND Meshgrid techniques')
    (data_meshed, param_valid) = space_parameters_ND(eTotal, phiEff, eMax, eMin, eStep,
                                                                     phiMax, phiMin, phiStep,
                                                                     eMaxLast, eMinLast,
                                                                     phiMaxLast, phiMinLast)
    #(E0, PHI0, E1, PHI1, E2, PHI2, E3, PHI3, ...) = data_meshed
    data_valid = tuple(item[param_valid] for item in data_meshed)
    for item in zip(*data_valid):
        print(item)
    param_count = np.count_nonzero(param_valid)
    print('Number of combinations:                            ', param_count)
    NonZeroValues = param_valid
    for item in data_meshed:
        NonZeroValues = np.logical_and(NonZeroValues, item != 0.0)
    paramNonZero_count = np.count_nonzero(NonZeroValues)
    print('Number of combinations with all non-zero elements: ', paramNonZero_count)
