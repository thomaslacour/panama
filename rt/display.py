import numpy as np
import matplotlib.pyplot as plt
import math as m

from . import waves as waves
from rt import coefficients as compute

from material import db as mat


def init_lp(data, **kwargs):
    """
    Initialization of layer properties.
    lp = (d, CL, CT, rhoL, rhoT, alphaL, alphaT)
    """
    d = []
    CL = []; CT = []
    rhoL = []; rhoT = []
    alphaL = []; alphaT = []
    omega = kwargs.get('omega', 0)
    # layers
    for layer in data['layers'].values():
        d.append(layer[1])
        mat_name = layer[0]
        if mat_name=='meta':
            metaparam = layer[2:]
            prop = mat.properties(mat_name, omega/2/np.pi, *metaparam)
            # prop = [u[0] for u in prop] # avoid array (1,)
        else:
            prop = mat.properties(mat_name, omega/2/np.pi)
        # prop = [ u[0] for u in prop if u!=[] ] # avoid array (1,)
        if prop==[]:
            return
        CL.append(prop[0])
        CT.append(prop[1])
        rhoL.append(prop[2])
        rhoT.append(prop[3])
        alphaL.append(prop[4])
        alphaT.append(prop[5])
    lp = (d, CL, CT, rhoL, rhoT, alphaL, alphaT)
    mat_name = data['halfspace_left']
    prop = mat.properties(mat_name, omega/2/np.pi)
    # prop = [u[0] for u in prop] # avoid array (1,)
    [ lp[i+1].insert(0, prop[i]) for i in range(6) ]
    mat_name = data['halfspace_right']
    prop = mat.properties(mat_name, omega/2/np.pi)
    # prop = [u[0] for u in prop] # avoid array (1,)
    [ lp[i+1].append(prop[i]) for i in range(6) ]
    return lp


def matrix_dim(CT):
    """
    DIMENSION of THE SCATTERING MATRIX
    test for left and right halfspace CT
    """
    if CT[0]==0 and CT[-1]==0:
        return 2
    elif CT[0]!=0 and CT[-1]!=0:
        return 4


def transfert_matrix_versus_omega(data): # OP
    f = np.linspace(data['f_min'],data['f_max'],data['f_num'])
    f = f[ f!=0 ]
    lp = init_lp(data, omega=2*np.pi*f[0])
    N = matrix_dim(CT=lp[2])
    M = np.ndarray((N,N,f.size),complex)
    name_half_space = data['halfspace_left']
    c = lp[1][0]
    theta = m.radians(data['theta_fix'])
    s1 = m.sin(theta)/c
    for i, omega in enumerate(2.0*np.pi*f):
        # lp = init_layers_prop(data, lm, layers_prop, omega=omega)
        lp = init_lp(data, omega=omega)
        M[:,:,i] = waves.transfert_matrix(omega, s1, data['unit_cells'], *lp)
    fig, axs = plt.subplots(ncols=N, nrows=N)
    for i in range(N):
        for j in range(N):
            axs[i,j].plot(f,np.real(M[i,j,:]),color='blue')
            axs[i,j].plot(f,np.imag(M[i,j,:]),color='red')
    plt.tight_layout()


def scattering_matrix_versus_omega(data): # OP
    f = np.linspace(data['f_min'],data['f_max'],data['f_num'])
    f = f[ f!=0 ]
    lp = init_lp(data, omega=2*np.pi*f[0])
    CL_left = lp[1][0]
    N = matrix_dim(CT=lp[2])
    CT = lp[2]
    if CT[0]==0 and CT[-1]==0:
        t = lambda i: 'L'
    elif CT[0]!=0 and CT[-1]!=0:
        t = lambda i: 'T' if i%2 else 'L'
    s = lambda i: '+' if i<N/2 else '-'
    S = np.ndarray((N,N,f.size),complex)
    theta = m.radians(data['theta_fix'])
    s1 = m.sin(theta)/CL_left
    for i, omega in enumerate(2.0*np.pi*f):
        lp = init_lp(data, omega=omega)
        S[:,:,i] = waves.scattering_matrix(omega, s1, data['unit_cells'], *lp)
    fig1, ax1 = plt.subplots(ncols=N, nrows=N)
    fig1.suptitle('Scattering matrix |Sij|', fontsize=14, fontweight='bold')
    fig2, ax2 = plt.subplots(ncols=N, nrows=N)
    fig2.suptitle('Scattering matrix arg(Sij)', fontsize=14, fontweight='bold')
    for i in range(N):
        for j in range(N):
            ax1[i,j].plot(f,np.abs(S[i,j,:]),color='black')
            ax1[i,j].set_title('S'+str(i//2+1)+str(j//2+1)+' '+t(j)+s(j)+' -> '+t(i)+s(i),fontsize=10, fontweight='bold')
            ax1[i,j].tick_params(axis='both', labelsize=8) 
            ax2[i,j].plot(f,np.angle(S[i,j,:])/np.pi,color='black')
            ax2[i,j].set_title('S'+str(i//2+1)+str(j//2+1)+' '+t(j)+s(j)+' -> '+t(i)+s(i),fontsize=10, fontweight='bold')
            ax2[i,j].tick_params(axis='both', labelsize=8)
            if j==0:
                ax1[i,j].set_ylabel('modulus', fontsize=8)
                ax2[i,j].set_ylabel('phase/$\pi$', fontsize=8)
            if i==3:
                ax1[i,j].set_xlabel('frequency (MHz)', fontsize=8)
                ax2[i,j].set_xlabel('frequency (MHz)', fontsize=8)

    pad, w_pad, h_pad = 1.5, 0.05, 0.05
    fig1.tight_layout(pad=pad, w_pad=w_pad, h_pad=h_pad)
    fig2.tight_layout(pad=pad, w_pad=w_pad, h_pad=h_pad)



def rt_versus_omega(data):
    freq_vec     = np.linspace(data['f_min'], data['f_max'], data['f_num'])
    theta        = [ data['theta_fix'] ]
    left         = (data['halfspace_left'], None) # None thickness
    right        = (data['halfspace_right'], None) # None thickness
    lay          = list(data['layers'].values()) # unpacked layers configuration
    stack_config = [left, *lay, right]
    stack_list   = compute.stacks(freq_vec, stack_config) # make a stack
    R, T = compute.rt(freq_vec, theta, stack_list , pola='LL')
    fig, ax = plot1d(freq_vec, R, T, phase=False)
    plt.xlabel('frequency (MHz)', fontsize=12)
    title = get_infos_for_title(data)
    fig.suptitle(title, fontsize=14, fontweight='bold')

    return freq_vec, R, T



def rt_versus_theta(data):
    theta_vec    = np.linspace(data['theta_min'], data['theta_max'], data['theta_num'])
    freq         = [ data['f_fix'] ]
    left         = (data['halfspace_left'], None) # None thickness
    right        = (data['halfspace_right'], None) # None thickness
    lay          = list(data['layers'].values()) # unpacked layers configuration
    stack_config = [left, *lay, right]
    stack_list   = compute.stacks(freq, stack_config)
    R, T = compute.rt(freq, theta_vec, stack_list , pola='LL')
    fig, ax = plot1d(theta_vec, R, T, phase=False)
    plt.xlabel('Angle (°C)', fontsize=12)
    title = get_infos_for_title(data)
    fig.suptitle(title, fontsize=14, fontweight='bold')

    return theta_vec, R, T



def get_infos_for_title(data):
    if len(data['layers'].keys()) == 1:
        key = [ u for u in data['layers'] ][0]
        layerprop = data['layers'][key]
        if layerprop[0]=='meta':
            layername = data['layers'][key][0] + '('
            layername = layername \
                    + str(layerprop[1]) + ' mm,' \
                    + str(layerprop[2]) + '+' + str(layerprop[3]) + ', ' \
                    + str(layerprop[4]) + ", " \
                    + str(layerprop[5]) + ", " \
                    + str(layerprop[6]) + ")"
            layername = '| ' + layername + ' |'
        else:
            layername = data['layers'][key][0] + '(' + str(layerprop[1]) + ' mm)'
    else:
        layername = "░░"

    return '{} {} {}' . format(data['halfspace_left'],layername,data['halfspace_right'])



def plot1d(x, R, T, phase=False):
    if phase==False:
        fig1, ax1 = plt.subplots(ncols=1, nrows=1)
        ax1 = np.array([ax1])
    else:
        fig1, ax1 = plt.subplots(ncols=1, nrows=2)

    # module :::
    ax1[0].plot(x,R,color='magenta', label=r'$R_{LL}$')
    ax1[0].plot(x,T,color='cyan',label=r'$T_{LL}$')
    ax1[0].set_ylabel('modulus', fontsize=12)
    ax1[0].plot\
    (
        x,np.ones(R.shape)-T-R,
        color='gray',
        label=r'$1-T_{LL}-R_{LL}$',
    )
    ax1[0].legend()

    # phase :::
    if phase==True:
        ax1[1].plot(x,np.angle(R)/np.pi,color='magenta')
        ax1[1].plot(x,np.angle(T)/np.pi,color='cyan')
        ax1[1].set_ylabel('phase/$\pi$', fontsize=12)

    plt.ylim(0,1)
    plt.yscale("linear")
    return fig1, ax1



def plot2d(x, y, R, T, phase=False):
    if phase==True:
        fig1, ax = plt.subplots(ncols=2, nrows=2)
        xylim = [x[0], x[-1], y[0], y[-1]]
        a = x[-1]/y[-1]
        cm = 'gnuplot2'
        cm2 = 'twilight'
        ax[0,0].imshow(abs(R)**2, extent=xylim, aspect=a, origin = 'lower', cmap=cm)
        ax[0,0].set_title('|$R_{LL}$|²')
        ax[1,0].imshow(np.angle(R)/np.pi, extent=xylim, aspect=a, origin = 'lower', cmap=cm2)
        ax[1,0].set_title('arg($R_{LL}$)')
        ax[0,1].imshow(abs(T)**2, extent=xylim, aspect=a, origin='lower', cmap=cm)
        ax[0,1].set_title('|$T_{LL}$|²')
        ax[1,1].imshow(np.angle(T)/np.pi, extent=xylim, aspect=a, origin = 'lower', cmap=cm2)
        ax[1,1].set_title('arg($T_{LL}$)')
        for i in range(2): ax[1,i].set_xlabel('Angle (°)', fontsize=12)
        for i in range(2): ax[i,0].set_ylabel('Frequency (MHz)', fontsize=12)

    elif phase==False:
        fig1, ax = plt.subplots(ncols=2, nrows=1)
        xylim = [x[0], x[-1], y[0], y[-1]]
        a = x[-1]/y[-1]
        cm = 'gnuplot2'
        cm2 = 'twilight'
        Rplot = ax[0].imshow(abs(R)**2, extent=xylim, aspect=a, origin = 'lower', cmap=cm)
        ax[0].set_title('|$R_{LL}$|²', y=1.3)
        ax[0].set_ylabel('Frequency (MHz)', fontsize=12)
        # ax[1].imshow(abs(T)**2, extent=xylim, aspect=a, origin='lower', cmap=cm)
        # ax[1].set_title('|$T_{LL}$|²')
        Aplot = ax[1].imshow(1-abs(T)**2-abs(R)**2, extent=xylim, aspect=a, origin='lower', cmap=cm)
        ax[1].set_yticklabels(['']*10)
        ax[1].set_title('1-|$R_{LL}$|²-|$T_{LL}$|²', y=1.3)
        for i in range(len(ax)): ax[i].set_xlabel('Angle (°)', fontsize=12)
        cax = fig1.add_axes([.125, .75, .35, .03])
        fig1.colorbar(Rplot, ax=ax[0], orientation='horizontal', fraction=0.044, cax=cax)
        cax.xaxis.tick_top()
        cax = fig1.add_axes([.547, .75, .35, .03])
        fig1.colorbar(Aplot, ax=ax[1], orientation='horizontal', fraction=0.044, cax=cax)
        cax.xaxis.tick_top()
    return ax






def PlotRT(CSVfile, title=""):
  # Trace les coefficients R et T depuis un fichier csv

  # lecture et extraction des données
  # data = readtable(CSVfile);
  # labels = data.Properties.VariableNames;
  # x = data{:,1};
  # R = data{:,2};
  # T = data{:,3};
  # A = sqrt(1-R.^2-T.^2);

  # % affichage
  # FigCoeff = figure;

  # plot(x, R.^2, 'linewidth', 2, 'color', [0.9,0.3,0.5]), hold on
  # plot(x, T.^2, 'linewidth', 2, 'color', [0.3,0.7,0.8])

  # plot(x, A.^2, 'linewidth', 2, 'color', 0.5*[1,1,1]), hold off

  # set(gca, 'YScale', 'log')

  # % annotations
  # if strcmp(labels{1},'freq_MHz')
  #   x = x*1e-3;
  #   xlabel('Fréquence (kHz)')
  # elseif strcmp(labels{1},'angle_deg')
  #   xlabel('Angle (degré)')
  # end
  # ylabel('Module')
  # legend('|R|²','|T|²', '1-|R|²-|T|²')
  # plt.title(varargin{1})
  pass
