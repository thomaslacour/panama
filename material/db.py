from .mat import *

def properties(Mat_n, freq, *args):
    """
    MATERIAL PROPERTIES
    Mat_n: name of the material
    freq: frequency (scalar or vector)
    *args: (mat, inc, r, phi, poly). See HOMOGENEOUSPROP().

    return phase velocity, attenuation coefficient and density
    """

    freq = freqformat(freq)


    # heterogeneous material.
    # ------------------------------------
    #   special case that call the program mst.py
    #   to compute effective properties for an
    # ------------------------------------
    if Mat_n in ['meta']:
        mat, inc = args[0], args[1]
        r = args[2]
        phi = args[3]
        poly = args[4]
        try:    pola = args[5]
        except: pola = None
        (cL, cT, rhoL, rhoT, alphaL, alphaT) = homogeneousprop(freq, mat, inc, r=r, phi=phi, poly=poly, pola=pola)

    # external material.
    # ------------------------------------
    #   an external program can use regex pattern 
    #   to replace values and temporarily add a material
    # ------------------------------------
    elif Mat_n in [ None ]:
        rhoL    =  None
        cL      =  None
        alphaL  =  None
        cT      =  None
        alphaT  =  None

    # ------------------------------------
    # the folowing material are homogenous
    # ------------------------------------
    elif Mat_n in ['eau', 'water', 'h2o']:
        rhoL    =  1.00
        cL      =  1.5
        alphaL  =  0.0
        cT      =  0.0
        alphaT  =  0.0
    elif Mat_n in ['air']:
        ### cf Onda_Gases ("Air - at 20°C")
        rhoL    =  0.001293
        cL      =  0.344
        alphaL  =  0.0
        cT      =  0.0
        alphaT  =  0.0
    elif Mat_n in ['PU', 'panama_PU']:
        rhoL    =  1.02
        # alphaL  =  0.001*freq**(1.4)
        # cL      =  1.52
        # alphaT  =  200 * alphaL
        # cT      =  0.18
        alphaL  =  0.05*freq**(1) # mesure par acoustique
        cL      =  1.5 # mesure par acoustique
        alphaT  =  6.3*freq**0.8 # mesure par DMA (naval group)
        cT      =  0.072*np.exp((freq)**0.2) # mesure par DMA (naval group)
    elif Mat_n in ['acier', 'steel']:
        ### cf Onda_Solids ("Steel - mild")
        rhoL    =  7.8
        cL      =  5.9
        alphaL  =  0.0 * freq
        cT      =  3.2
        alphaT  =  0.0 * freq
    elif Mat_n in ['aluminium', 'alu', 'aluminum']:
        ### cf Onda_Solids ("Steel - mild")
        rhoL    =  2.7
        cL      =  6.32
        alphaL  =  0.0 * freq
        cT      =  3.13
        alphaT  =  0.0 * freq
    elif Mat_n in ['Silicone_poreux', 'PDMS_poreux', 'PDMS', 'billes', 'panama_billes']:
        rhoL    =  0.70 # valeur romain
        cL      =  0.2+0.1*freq # mesure acoustique
        alphaL  =  9*freq**(0.9) # mesure acoustiquee
        # ??? supposition !!!
        cT      =  0.05
        alphaT  =  18 * alphaL

    ###################################################################
    # LIST OF USER'S MATERIALS
    # ------------------------------------
    #   just copy pas the pattern of code.
    #   It needs at least the entries rhoL, cL, alphaL, cT and alphaT
    # elif Mat_n in [ 'material_name' ]:
    #     rhoL    = 0
    #     cL      = 0
    #     alphaL  = 0
    #     cT      = 0
    #     alphaT  = 0
    ###################################################################


    if not 'rhoT' in locals() :
        rhoT=rhoL

    return formatprop(Mat_n, freq, cL, cT, rhoL, rhoT, alphaL, alphaT)
