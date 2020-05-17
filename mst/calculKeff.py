import numpy as np
from numpy import imag, real, exp # -------------------- complexe ----- |
from numpy import pi, sqrt # --------------------------- math --------- |
# ------------------------------------------------------ array -------- |
from numpy import ndarray, array, arange, zeros, gradient, linspace
import matplotlib.pyplot as plt # ---------------------- figure ------- |

from .calculScatteringCoeffsSn import *
from material import db as mat

def calculKeff(freq, Model, nMode, Mat, Inc, r_mean_mm, poly, phivol, **kwargs):
    """
        Model       = choix du modele MST (Foldy, Waterman/Truell ou Linton/Martin)
        nMode       = nbre de modes à considerer pour le calcul de fonct. de  diff.
                      par exemple, n=2 ---> monopolaire + dipolaire
        Mat         = choix du materiau pour la matrice
        Inc         = choix du materiau pour les inclusions
        r_mean_mm   = rayon moyen des inclusions en mm
        poly        = coefficient de variation de la distribution en taille
        phivol      = fraction volumique des inclusions en %
        d           = epaisseur de l'echantillon en mm
    """

    verbose = kwargs.get('verbose', 1)

    if not 'polar' in kwargs:
        print('No polarisation has been given. Longitudinal will be used as default')
        pola = 'L'
    else:
        pola = kwargs.get('polar')

    if pola == 'L':
        calculScatteringCoeffsSn = SnLL
    elif pola == 'T':
        calculScatteringCoeffsSn = SnTT

    freq = mat.freqformat(freq)

    # compute material properties
    rho0, lambda0, mu0, rho1, lambda1, mu1 = mat.initprop(freq,Mat,Inc,pola,**kwargs)

    poly = poly/100
    r_sigma_mm = r_mean_mm*poly
    NFREQ = len(freq)
    phivol = phivol/100
    omega = 2*pi*freq

    # elastic modulus
    if pola == 'L':
        M0 = (lambda0 + 2*mu0)
    if pola == 'T':
        M0 = mu0

    # wave number
    k0 = omega*sqrt(rho0/M0)
    # phase velocity
    c0 = omega/real(k0)
    # attenuation
    alpha0 = imag(k0)
    # indice
    n0 = k0/k0
    # impedance
    Z0 = rho0*omega/k0

    # computation fo all thz size
    r_min_1 = r_mean_mm - 3*r_sigma_mm
    r_max_1 = r_mean_mm + 3*r_sigma_mm
    Nb_r = 100 #nombre de rayons pour la distribution en taille
    Da = (r_max_1-r_min_1)/Nb_r
    if poly>0:
        radiusvaluesArray = arange(r_min_1, r_max_1, Da)
    elif poly==0:
        radiusvaluesArray = array([r_mean_mm])

    # initialization of scattering functions 0 and pi
    nradius = len(radiusvaluesArray)
    f_shape = (NFREQ,nradius)

    f0_f_radius = zeros( f_shape, dtype=complex )
    fpi_f_radius = zeros( f_shape, dtype=complex)

    nmax = nMode-1

    Sn_f_a_n = zeros( (NFREQ, nradius, nmax+1), dtype=complex )

    if verbose==1 and nradius > 1:
        import misc.progressbar as bar
        print('Computing diffusion functions from {} to {} µm: ' . format(r_min_1*1000,r_max_1*1000))
        pref='Progress'
        bar.printProgressBar(0, nradius, prefix = pref, suffix = '', length = 30)

    for itaille in range(nradius):
        rd = radiusvaluesArray[itaille]
        S, nmax = calculScatteringCoeffsSn(lambda0,mu0,rho0,lambda1,mu1,rho1,freq,rd,nmax)
        # Sn_f_a_n[:,itaille-1,:] = S #<----unused
        if pola=='L':
            TnLL = S[0,:,:]
            f0_f_radius[:,itaille-1] = calcul_f_sc(0, freq, k0, TnLL, nmax).T
            fpi_f_radius[:,itaille-1] = calcul_f_sc(pi, freq, k0, TnLL, nmax).T

        elif pola=='T':
            TnTT = S[0,:,:]
            tnTT = S[4,:,:]
            f0_f_radius[:,itaille-1] = f_tt(0, freq, k0, [TnTT, tnTT], nmax)
            fpi_f_radius[:,itaille-1] = f_tt(pi, freq, k0, [TnTT, tnTT], nmax)

        if verbose==1 and nradius > 1:
            # print(rd*1e3) # -------
            bar.printProgressBar(itaille + 1, nradius, prefix = pref, suffix = '', length = 30)

    radius_range = np.logical_and(radiusvaluesArray >= r_min_1, radiusvaluesArray <= r_max_1)

    f0_f_rd = f0_f_radius[:, radius_range]
    fpi_f_rd = fpi_f_radius[:, radius_range]

    radiusDistribArray = radiusvaluesArray[radius_range]

    # distribution gaussienne en n(a)
    if poly > 0:
        GaussianArray = exp(-0.5*((radiusDistribArray-r_mean_mm)/r_sigma_mm)**2)
        GaussianArrayNorm = (1/sum(GaussianArray))*GaussianArray
        N_V = phivol / sum( GaussianArrayNorm * ( (4/3)*pi*(radiusDistribArray**3) ) )
        na = N_V*GaussianArrayNorm
        NSCATTERERS_VOLDistrib = array([na]).T

    # monodisperse
    elif poly == 0:
        dirac = zeros( radiusDistribArray.shape )
        dirac[(radiusvaluesArray == r_mean_mm)] = 1
        GaussianArrayNorm = dirac
        N_V = phivol/sum( GaussianArrayNorm*( (4/3)*pi*(radiusDistribArray**3) ) )
        na = N_V*GaussianArrayNorm
        NSCATTERERS_VOLDistrib = array([na]).T

    model = Model.lower()
    # Foldy
    #   keff^2 = k0^2 + 4*pi*n(a)*f(0)
    if model in ['foldy']:
        keff2 = k0**2 + 4*pi*np.einsum('ji,i->j', f0_f_rd, na)

    # Watterman and Truell
    #   keff^2 = [k0 + 2*pi*n(a)*f(0)/k0]^2 - [2*pi*n(a)*f(pi)/k0]^2
    elif model in ['wt']:
        wt1 = 1 + (2*pi/k0**2)*np.einsum('ji,i->j', f0_f_rd - fpi_f_rd, na)
        wt2 = 1 + (2*pi/k0**2)*np.einsum('ji,i->j', f0_f_rd + fpi_f_rd, na)
        keff2 = k0**2*wt1*wt2

    # Linton and Martin
    # + Delta2_LM_polydisp(k0,NSCATTERERS_VOLDistrib, Sn_f_a_n,freq) ;
    # TODO: redo the matlab code

    k1 = real(sqrt(keff2))
    k2 = imag(sqrt(keff2))

    # we impose a positive attenuation and thus we fix the direction for the
    # propagation of the energy
    k2[k1<0] *= -1
    k1[k1<0] *= -1

    keff = k1 + 1j*k2
    neff = keff/k0
    ceff = omega/real(keff)
    Alphaeff = imag(keff)

    # group velocity
    # cgr = gradient(omega,real(keff)) ;
    # ou cgr = c0 ./ ( real(neff) + omega.*gradient(real(neff),omega)) ;

    # C. Aristegui (Wave Motion 47 (2010) 199-204)
    # rhoeff  = rho0*( 2*pi/k0^2 )*( (f0_f_rd - fpi_f_rd) * n(a) )
    rhoeff = np.einsum('ji,i->j', f0_f_rd - fpi_f_rd, na)
    rhoeff = np.einsum('i,i->i', 2*pi/k0**2, rhoeff)
    rhoeff+= rho0

    # Meff  = M0 / ( 1 + (2*pi/k0^2)*( (f0_f_rd + fpi_f_rd) * n(a) ) )
    Meff = np.einsum('ji,i->j', f0_f_rd + fpi_f_rd, na)
    Meff = np.einsum('i,i->i', 2*pi/k0**2, Meff)
    Meff = M0/(1+Meff)

    Zeff = rhoeff * omega / keff

    out = dict()
    out= {  'freq'          : freq,
            'k_eff'         : keff,
            'alpha_eff'     : Alphaeff,
            'c_eff'         : ceff,
            'n_eff'         : neff,
            'rho_eff'       : rhoeff,
            'M_eff'         : Meff,
            'Z_eff'         : Zeff,
            'radii'         : radiusDistribArray,
            'NScatVolDis'   : NSCATTERERS_VOLDistrib,
       }

    return out


# end of calculKeff.py ----
