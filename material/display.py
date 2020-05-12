import numpy as np
import math as m
import matplotlib.pyplot as plt

import matplotlib.style
import matplotlib as mpl
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['font.size'] = 12
mpl.rcParams['figure.titlesize'] = 'medium'
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'

from . import db as mat

def prop(mp, freq, **kwargs):

    mat_name = mp[0]
    d = mp[1]
    freq = freq[freq!=0]
    if mat_name == 'meta':
        (cL, cT, rhoL, rhoT, alphaL, alphaT) = mat.homogeneousprop(freq, *mp[2:], **kwargs)
        matprop = mat.properties(mp[2], freq)
    else:
        (cL, cT, rhoL, rhoT, alphaL, alphaT) = mat.properties(mat_name, freq)
        pass

    mu_, lambda_ =  mat.lamecoeff(freq, cL, cT, alphaL, alphaT, rhoL, rhoT)

    K = lambda_ + 2/3*mu_
    ML = K + 4/3*mu_
    MT = mu_

    omega = 2*np.pi*freq
    kL=omega/cL + 1j*alphaL
    kT=omega/cT + 1j*alphaT
    ZL=rhoL*omega/kL
    ZT=rhoT*omega/kT

    #init display
    fig, ax_real = plt.subplots(ncols=4, nrows=2, sharex=True)

    color_real = 'black'
    color_imag = 'gray'
    linestyle_real = 'solid'
    linestyle_imag = 'dashed'

    ax_imag = ax_real.copy()
    for i in range(ax_imag.shape[0]):
        for j in range(ax_imag.shape[1]):
            ax_imag[i,j] = ax_real[i,j].twinx()
            ax_imag[i,j].tick_params(axis='y', labelcolor=color_imag)

    # cL_ax     = plt.subplot(grid[0,0])
    # cT_ax     = plt.subplot(grid[0,1])
    # alphaL_ax = plt.subplot(grid[1,0])
    # alphaL_ax = cL_ax.twinx()
    # alphaT_ax = cT_ax.twinx()
    # ZL_ax     = plt.subplot(grid[2,0])
    # ZT_ax     = plt.subplot(grid[2,1])
    # rho_ax    = plt.subplot(grid[0,2])
    # M_ax      = plt.subplot(grid[1,2])
    # mu_ax     = plt.subplot(grid[2,2])

    x = freq

    fig.subplots_adjust(wspace=0.2,hspace=0.15)

    plt.axis('off')

    # wavenumber
    ax_real[0,0].plot(x, cL, color=color_real, linestyle=linestyle_real)
    ax_imag[0,0].plot(x, alphaL, color=color_imag, linestyle=linestyle_imag)
    ax_real[1,0].plot(x, cT, color=color_real, linestyle=linestyle_real)
    ax_imag[1,0].plot(x, alphaT, color=color_imag, linestyle=linestyle_imag)
    ax_real[0,0].set_title(r'$k=\omega/c+i\alpha$')
    # ax.tick_params(axis="x",direction="in", pad=-15)
    # density
    ax_real[0,1].plot(x, np.real(rhoL), color=color_real, linestyle=linestyle_real)
    ax_imag[0,1].plot(x, np.imag(rhoL), color=color_imag, linestyle=linestyle_imag)
    ax_real[1,1].plot(x, np.real(rhoT), color=color_real, linestyle=linestyle_real)
    ax_imag[1,1].plot(x, np.imag(rhoT), color=color_imag, linestyle=linestyle_imag)
    ax_real[0,1].set_title(r'$\rho$')
    # impedence
    ax_real[0,2].plot(x, np.real(ZL), color=color_real, linestyle=linestyle_real)
    ax_imag[0,2].plot(x, np.imag(ZL), color=color_imag, linestyle=linestyle_imag)
    ax_real[1,2].plot(x, np.real(ZT), color=color_real, linestyle=linestyle_real)
    ax_imag[1,2].plot(x, np.imag(ZT), color=color_imag, linestyle=linestyle_imag)
    ax_real[0,2].set_title(r'$Z=\rho\omega/k$')
    # elastic modulus
    ax_real[0,3].plot(x, np.real(ML), color=color_real, linestyle=linestyle_real)
    ax_imag[0,3].plot(x, np.imag(ML), color=color_imag, linestyle=linestyle_imag)
    ax_real[1,3].plot(x, np.real(MT), color=color_real, linestyle=linestyle_real)
    ax_imag[1,3].plot(x, np.imag(MT), color=color_imag, linestyle=linestyle_imag)
    ax_real[0,3].set_title(r'$M_L = K+\frac{4}{3}\mu$')
    ax_real[1,3].set_title(r'$M_T=\mu$')
    # axis
    ax_real[0,0].set_ylabel(r'Longitudinal')
    ax_real[1,0].set_ylabel(r'Transversal')
    # [ ax_real[1,i].set_xlabel(r'Frequency (MHz)') for i in range(ax_real.shape[1]) ]

    ticks = kwargs.get('ticks', True)
    for a in [ax_real, ax_imag]:
        for axi in a.flat:
            if not ticks:
                axi.set_xticklabels([])
                axi.set_yticklabels([])
            axi.xaxis.set_major_locator(plt.MaxNLocator(3))
            axi.yaxis.set_major_locator(plt.MaxNLocator(3))
            # axi.ticklabel_format(axis="y",style='sci',scilimits=(1,4))
            # axi.tick_params(axis="y",direction="in", pad=-25)
            axi.autoscale(enable=True, axis='both', tight=True)
            axi.tick_params(axis='both', which='major', labelsize=8)

    fig.text(0.5, 0.04, 'Frequency (MHz)', ha='center')

    # fig.tight_layout()


    return (cL, alphaL, cT, alphaT)
