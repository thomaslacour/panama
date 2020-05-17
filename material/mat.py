import numpy as np
import material.db as db

def freqformat(freq):
    """
    Manage the right format for frequency input (scalar or vector)
    return an array every time of size (1,) or (N,).
    """

    if type(freq) is list:
        freq=np.array(freq)
    elif type(freq) in  [int, float, np.float64]:
        # freq = np.array( (freq,) )
        # freq = freq
        pass

    if type(freq) is np.ndarray:
        freq[freq==0]=1e-9

    return freq

def initprop(freq, Mat_n, Inc_n, pola, **kwargs):
    """
    Init material properties in term of rho (density), lambda (first 
    lamé coefficient) and mu (second lamé coefficient).
    Also return the right format for "freq"
    """

#     # # avoid zero frequency for computation
#     freq = freq[ freq!=0 ]

    # material properties
    cL0, cT0, rhoL0, rhoT0, alphaL0, alphaT0 = db.properties(Mat_n, freq)
    cL1, cT1, rhoL1, rhoT1, alphaL1, alphaT1 = db.properties(Inc_n, freq)

    # switching to complex values of mu & lambda
    mu0, lambda0 =  lamecoeff(freq, cL0, cT0, alphaL0, alphaT0, rhoL0, rhoT0)
    mu1, lambda1 =  lamecoeff(freq, cL1, cT1, alphaL1, alphaT1, rhoL1, rhoT1)

    if   pola=='L': rho0,rho1=rhoL0,rhoL1
    elif pola=='T': rho0,rho1=rhoT0,rhoT1

    return rho0, lambda0, mu0, rho1, lambda1, mu1



def missing_material_warning(Mat_n):
    print('WARNING : Material named \'{}\' not found in db.'.format(Mat_n))
    try:
        import sys, re, os
        path = os.getcwd() + '/material/db.py'
        print('Please add it to the database {}, or choose it within the following list :\n' . format(path))
        with open(path, 'r') as file: data = file.read()
        [ print(i) for i in re.findall(r'[M]at_n in (.*):', data)[:-1] ]
    except:
        pass
    finally:
        sys.exit(0)
    return



def formatprop(Mat_n, freq, cL=None, cT=None, rhoL=None, rhoT=None, alphaL=None, alphaT=None):

    u = freq
    if type(u) == list and len(u) == 1: u = 1
    if type(u) == np.ndarray: u = np.ones( u.shape )

    try:
        return (cL*u, cT*u, rhoL*u, rhoT*u, alphaL*u, alphaT*u)
        # if not 'rhoT' in locals() :
        #     # return (cL, cT, rhoL, rhoL, alphaL, alphaT)
        #     return (cL*u, cT*u, rhoL*u, rhoL*u, alphaL*u, alphaT*u)
        # else:
        #     # return (cL, cT, rhoL, rhoT, alphaL, alphaT)
        #     return (cL*u, cT*u, rhoL*u, rhoT*u, alphaL*u, alphaT*u)
    except:
        print('WARNING : Material named \'{}\' not found in db.'.format(Mat_n))
        try:
            import sys, re, os
            path = os.getcwd() + '/material/db.py'
            print('Please add it to the database {}, or choose it within the following list :\n' . format(path))
            with open(path, 'r') as file: data = file.read()
            [ print(i) for i in re.findall(r'Mat_n in (.*):', data) ]
        except:
            pass
        finally:
            sys.exit(0)


def lamecoeff(freq, cL, cT, alphaL, alphaT, rhoL, rhoT=None):
    """
    Lamé coefficient
    """
    kL = wavenumber(freq, cL, alphaL)
    kT = wavenumber(freq, cT, alphaT)
    mu0 = rhoT*(freq*2*np.pi)**2 / kT**2
    lambda0 = rhoL*(freq*2*np.pi)**2 / kL**2 - 2*mu0
    return mu0, lambda0


def wavenumber(freq, c, alpha):
    """
    Use phase velocity and attenuation to return the complex
    wave number.
    """
    freq = freqformat(freq)
    k = freq*2*np.pi/c + 1j*alpha
    if np.imag(np.array(k)).all()==0:
        k = np.real(k)
    return k


def homogeneousprop(freq, mat, inc, r, phi, poly, pola=None, **kwargs):
    """
    HOMOGENEOUS PROPERTIES using multiple scattering theory.
    freq: frequency
    mat: material name for the matrix
    inc: material name for the inclusion
    r: mean radius of inclusion (mm)
    phi: volume fraction of inclusion in %
    poly: polydispersity of the inclusions radii
    """
    freq = freqformat(freq)
    c, rho, alpha = [ {} for i in range(3) ]

    from mst.calculKeff import calculKeff
    # default values
    N = kwargs.get('nMode', 2)
    Model = kwargs.get('Model', 'WT')

    if pola==None:
        pola = ['L', 'T']

    for p in pola:
        out = calculKeff(
            Model     = Model,
            nMode     = N,
            fMAX      = 0,
            Mat       = mat,
            Inc       = inc,
            r_mean_mm = r,
            phivol    = phi,
            poly      = poly,
            d         = 0,
            freq      = freq,
            polar     = p,
            verbose   = 0,
            )
        out = (out['c_eff'], out['rho_eff'], out['alpha_eff'])
        # out = mst(freq, mat, inc, r, phi, poly, p, **kwargs)
        if (1,) in [np.array([freq]).shape, freq.shape]:
            c[p] = out[0][0]
            rho[p] = out[1][0]
            alpha[p] = out[2][0]
        else:
            c[p] = out[0]
            rho[p] = out[1]
            alpha[p] = out[2]

    if len(pola)==2:
        return (c['L'], c['T'], rho['L'], rho['T'], alpha['L'], alpha['T'])
    elif len(pola)==1:
        if pola == 'L':
            return (c['L'], 0, rho['L'], 0, alpha['L'], 0)
        elif pola == 'T':
            return (0, c['T'], 0, rho['T'], 0, alpha['T'])

