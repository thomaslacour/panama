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

    # freq = freqformat(kwargs.get('f', 1))
    # avoid zero frequency for computation
    # freq = freq[ freq!=0 ]

    freq,rho0,lambda0,mu0,rho1,lambda1,mu1 = mat.initprop(freq,Mat,Inc,pola,**kwargs)

    poly = poly/100
    r_sigma_mm = r_mean_mm*poly
    NFREQ = len(freq)
    phivol = phivol/100
    omega = 2*pi*freq

    if pola == 'L':
        M0 = (lambda0 + 2*mu0)
    if pola == 'T':
        M0 = mu0
    k0 = omega*sqrt(rho0/M0)
    c0 = omega/real(k0)
    alpha0 = imag(k0)
    n0 = k0/k0
    Z0 = rho0*omega/k0

    # calcul coeffs Sn pour chaque taille
    r_min_1 = r_mean_mm - 3*r_sigma_mm
    r_max_1 = r_mean_mm + 3*r_sigma_mm
    Nb_r = 100 #nombre de rayons pour la distribution en taille
    Da = (r_max_1-r_min_1)/Nb_r
    if poly>0:
        radiusvaluesArray = arange(r_min_1, r_max_1, Da)
    elif poly==0:
        radiusvaluesArray = array([r_mean_mm])

    nradius = len(radiusvaluesArray)
    f0_f_radius = zeros( (NFREQ, nradius), dtype=complex )
    fpi_f_radius = zeros( (NFREQ, nradius), dtype=complex)

    nmax = nMode-1

    Sn_f_a_n = zeros( (NFREQ, nradius, nmax+1), dtype=complex )

    if verbose==1:
        print('Computing diffusion functions from {} to {} µm: ' . format(r_min_1*1000,r_max_1*1000))
        pref='Progress'
        printProgressBar(0, nradius, prefix = pref, suffix = '', length = 30)

    for itaille in range(nradius):
        rd = radiusvaluesArray[itaille]
        S, nmax = calculScatteringCoeffsSn(lambda0,mu0,rho0,lambda1,mu1,rho1,freq,rd,nmax)
        # Sn_f_a_n[:,itaille-1,:] = S #<----pas utilisé !
        if pola=='L':
            TnLL = S[0,:,:]
            f0_f_radius[:,itaille-1] = calcul_f_sc(0, freq, k0, TnLL, nmax).T
            fpi_f_radius[:,itaille-1] = calcul_f_sc(pi, freq, k0, TnLL, nmax).T

        elif pola=='T':
            TnTT = S[0,:,:]
            tnTT = S[4,:,:]
            f0_f_radius[:,itaille-1] = f_tt(0, freq, k0, [TnTT, tnTT], nmax)
            fpi_f_radius[:,itaille-1] = f_tt(pi, freq, k0, [TnTT, tnTT], nmax)

        if verbose==1:
            # print(rd*1e3) # -------
            printProgressBar(itaille + 1, nradius, prefix = pref, suffix = '', length = 30)

    radius_range = (radiusvaluesArray >= r_min_1) & (radiusvaluesArray <= r_max_1)
    f0_f_rd = f0_f_radius[:, radius_range]
    fpi_f_rd = fpi_f_radius[:, radius_range]

    radiusDistribArray = radiusvaluesArray[radius_range]

    # distribution gaussienne en n(a)
    if poly>0:
        GaussianArray = exp(-0.5*((radiusDistribArray-r_mean_mm)/r_sigma_mm)**2)
        GaussianArrayNorm = (1/sum(GaussianArray))*GaussianArray
        N_V = phivol / sum( GaussianArrayNorm * ( (4/3)*pi*(radiusDistribArray**3) ) )
        NSCATTERERS_VOLDistrib = N_V*GaussianArrayNorm
        NSCATTERERS_VOLDistrib = array([NSCATTERERS_VOLDistrib]).T
    # monodisperse
    elif poly == 0:
        dirac = zeros( radiusDistribArray.shape )
        dirac[(radiusvaluesArray == r_mean_mm)] = 1
        GaussianArrayNorm = dirac
        N_V = phivol/sum( GaussianArrayNorm*( (4/3)*pi*(radiusDistribArray**3) ) )
        NSCATTERERS_VOLDistrib = N_V*GaussianArrayNorm
        NSCATTERERS_VOLDistrib = array([NSCATTERERS_VOLDistrib]).T

    # %% Calcul du vecteur d'onde effective
    #
    # switch Model
    #
    #     case 'Foldy'   % Modele de Foldy
    #         keff2 = k0.^2 + 4*pi     *     f0_f_rd          * NSCATTERERS_VOLDistrib ;
    #
    #     case 'WT'      % Modele de Waterman & Truell
    # %        keff2 = k0.^2 + 4*pi     *     f0_f_rd          * NSCATTERERS_VOLDistrib ...
    # %                      + (4*pi^2./k0.^2) .* ((f0_f_rd.^2-fpi_f_rd.^2) * NSCATTERERS_VOLDistrib.^2) ;
    #         keff2 = (k0.^2) .* ( ...
    #             (  ones(NFREQ,1) +    (2*pi./(k0.^2) ) .* ( (f0_f_rd - fpi_f_rd) *  NSCATTERERS_VOLDistrib)  )...
    #             .* (  ones(NFREQ,1) +    (2*pi./(k0.^2) ) .* ( (f0_f_rd + fpi_f_rd) *  NSCATTERERS_VOLDistrib)  )...
    #             );
    ####TODO: BAD SHAPE (N, 1) ? ############################# !
    k0 = k0.reshape( (len(k0),1) )
    omega = omega.reshape( (len(omega),1) )
    if type(rho0) is ndarray: rho0 = rho0.reshape( (len(rho0),1) )
    ########################################################## !
    keff2 = (k0**2) * (     \
        (  ones( (NFREQ,1) ) + (2*pi/(k0**2) ) * ( (f0_f_rd - fpi_f_rd) @ NSCATTERERS_VOLDistrib)  )         \
      * (  ones( (NFREQ,1) ) + (2*pi/(k0**2) ) * ( (f0_f_rd + fpi_f_rd) @ NSCATTERERS_VOLDistrib)  )      \
             )
    #
    #     case 'LM'      % Modele de Linton & Martin
    #         keff2 = k0.^2 + 4*pi     *     f0_f_rd          * NSCATTERERS_VOLDistrib ...
    #                       + Delta2_LM_polydisp(k0,NSCATTERERS_VOLDistrib, Sn_f_a_n,freq) ;
    #
    # end
    #
    #
    # % keff = sqrt(keff2) ;
    # k1 = real(sqrt(keff2)) ;
    k1 = real(sqrt(keff2))
    # k2 = imag(sqrt(keff2)) ;
    k2 = imag(sqrt(keff2))
    #
    # On impose l'atténuation positive ce qui revient à fixer la direction de
    # propagation de l'énergie pour un milieu passif
    for ii in range(keff2.shape[0]):
        if k2[ii] < 0:
            k1[ii] = - k1[ii]
            k2[ii] = - k2[ii]

    keff = k1 + 1j*k2
    neff = keff/k0
    ceff = omega/real(keff)
    Alphaeff = imag(keff)

    # Vitesse de groupe
    # cgr = gradient(omega,real(keff)) ;
    # ou cgr = c0 ./ ( real(neff) + omega.*gradient(real(neff),omega)) ;

    # C. Aristegui (Wave Motion 47 (2010) 199-204)
    # rhoeff  = rho0 .* ones(NFREQ,1) + ( 2*pi./ (k0.^2) )  .* ( (f0_f_rd - fpi_f_rd) * NSCATTERERS_VOLDistrib );
    rhoeff  = rho0*ones( (NFREQ,1) )+( 2*pi/(k0**2) ) * ( (f0_f_rd - fpi_f_rd) @ NSCATTERERS_VOLDistrib )
    # Meff  = M0 ./ ( ones(NFREQ,1) + ( 2*pi./ (k0.^2) )  .* ( (f0_f_rd + fpi_f_rd) * NSCATTERERS_VOLDistrib ) );
    Meff  = M0/( ones( (NFREQ,1) ) + ( 2*pi/(k0**2) ) * ( (f0_f_rd + fpi_f_rd) @ NSCATTERERS_VOLDistrib ) )
    Zeff = rhoeff * omega / keff
    Zw = 1.5*ones( Zeff.shape )

    # Transmission = (4*Zeff.*Zw./((Zw+Zeff).^2)) .* exp(i*keff*d);
    #### Transmission = (4*Zeff*Zw/((Zw+Zeff)**2)) * exp(1j*keff*d);
    # 
    # %%% estimation de keff en négligeant les interfaces du slab plongé dans l'eau
    # T = (4*Zeff.*Zw./((Zw+Zeff).^2)) ;
    # keffBis = log(T.*exp(i*keff*d))/i/d ;
    # ceffBis = omega./real(keffBis) ;
    # AlphaeffBis = imag(keffBis) ;
    # cgrBis = gradient(omega,real(keffBis)) ;
    #

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
    ##TODO: RESHAPE BAD RESHAPE !
    for u in out.keys():
        v = out[u]
        if len(v.shape)>1:
            if v.shape[1]==1:
                v = v.reshape( (v.shape[0],) )
        out[u] = v
    return out




# progress bar
import time
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()
# end of calculKeff.py ----
