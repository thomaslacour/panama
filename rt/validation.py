import numpy as np
import cmath as m
import matplotlib.pyplot as plt

from rt import actions as action
from rt import display as display

from material import db as material
from rt.coefficients import stacks, rt
from rt.display import plot1d

def test_1layer_NormalIncidence(*arg):
    print('\n1| Normal incidence for 1 layer immersed in a fluid')
    # --------------------------------------------------------
    # Test de validation du code en incidence normale.
    # --------------------------------------------------------
    lay = []
    # initialize the representation vectors
    freq_vec = np.linspace(0, 10, 301)
    theta = [0]
    # configuration of the layers, None thickness for halfspaces
    left  =   ('water', None)
    lay  += [ ('steel', 1) ]
    right =   ('water', None)
    stack_config = [left, *lay, right]
    # stacking the layers
    stack_list = stacks(freq_vec, stack_config)
    # compute RT
    R, T = rt(freq_vec, theta, stack_list , pola='LL')
    # display coefficients
    fig, ax = plot1d(freq_vec, R, T, phase=False)
    plt.xlabel('frequency (MHz)', fontsize=12)
    # save layers properties
    lp = [ tuple(u.tolist()) for u in stack_list[0] ]
    # save of numerical solution
    f_num, R_num, T_num = freq_vec, R, T #<--------
    # --------------------------------------------------------
    # formule analytique
    # --------------------------------------------------------
    # ARISTÉGUI, et al., Wave Motion (2010), Eq. (2)
    # --------------------------------------------------------
    h = lay[0][1]/2
    freqs = freq_vec
    rho1, cL1 = lp[3][1], lp[1][1]
    Z1 = rho1*cL1
    rho0, cL0 = lp[3][0], lp[1][0]
    Z0 = rho0*cL0
    k0 = 2*np.pi*freqs/cL0
    adZeff = Z1/Z0 # rapport d'impédance
    Keff = 2*np.pi*freqs/cL1
    Q = (1-adZeff)/(1+adZeff)
    buff = np.exp(4*1j*Keff*h)
    denominateur = 1-Q**2*buff
    T = (1-Q**2)*np.exp(2j*(Keff-k0)*h)/denominateur
    R = -Q*(1-buff)*np.exp(2j*k0*h)/denominateur
    # attention le centre du repère est au milieu de la plaque !
    a = np.exp(1j*k0*2*h) # <------!!
    R = R/a
    T = T*a
    f_ana, R_ana, T_ana = freqs, R, T#<----------

    plt.plot(freqs, np.abs(R)**2, linestyle='--', label='$|R_{LL}|^2~(th)$')
    plt.plot(freqs, np.abs(T)**2, linestyle='--', label='$|T_{LL}|^2~(th)$')
    # ax1[1].plot(freqs, np.angle(R)/np.pi, linestyle='--')
    # ax1[1].plot(freqs, np.angle(T)/np.pi, linestyle='--')
    plt.legend()

    # kd
    # ax2 = ax1[0].twiny()
    # f = np.array([0,0.5,1])*cL1/(2*h)
    # def tick_function(freq):
    #     x = freq*2*h/cL1
    #     return ["%.3f" % z for z in x]
    # ax2.set_xlim(ax1[0].get_xlim())
    # ax2.set_xticks(f)
    # ax2.set_xticklabels(tick_function(f))
    # ax2.set_xlabel("$d\,/\,\lambda$")
    # plt.show()

    return
    # return (f_num*2*h/cL1, R_num, T_num, f_ana*2*h/cL1, R_ana, T_ana)

def test_1layer_Angle(*arg):
    print('\n2| Incidence for 1 layer in a fluid')
    # --------------------------------------------------------
    # Test de validation du code sous incidence.
    # --------------------------------------------------------

    lay = []
    # initialize the representation vectors
    f, theta_vec = [1], np.linspace(0, 90, 901) # MHz, degrees
    # configuration of the layers, None thickness for halfspaces
    left  =   ('water', None)
    lay  += [ ('steel', 2) ]
    right =   ('water', None)
    stack_config = [left, *lay, right]
    # stacking the layers
    stack_list = stacks(f, stack_config)
    # compute RT
    R, T = rt(f, theta_vec, stack_list , pola='LL')
    # display coefficients
    fig, ax = plot1d(theta_vec, R, T, phase=False)
    plt.xlabel('angle (°)', fontsize=12)
    # save layers properties
    lp = [ tuple(u.tolist()) for u in stack_list[-1] ]
    # save of numerical solution
    theta_num, R_num, T_num = theta_vec, R, T #<--------

    #TODO: some bugs on analytical formula. The numerical result is good compare
    #      to previous test.
    return
    # formule analytique: 3(fluid) | 2(solid) | 1(fluid=3)
    # --------------------------------------------------------
    # BREKHOVSKIKH, Waves in Layered Media, p. 71
    # Eq. (10.7) et (10.8)
    # --------------------------------------------------------
    h = lp[0][1]
    f = f[0]
    # h = data['layers']['1'][1]
    theta = np.radians(theta_num)
    theta_num = np.degrees(theta.copy())
    c1 = lp[1][0]
    rho2 = lp[3][1]
    b2 = lp[2][1] # CT
    c2 = lp[1][1] # CL
    rho1 = lp[3][0]
    theta1 = theta.copy()
    arg = c2/c1*np.sin(theta1)
    arg[ arg>=1 ] = None
    theta2 = np.arcsin(arg)
    arg = b2/c1*np.sin(theta1)
    arg[ arg>=1 ] = None
    gamma2 = np.arcsin(arg)
    K2 = 2*np.pi*f/c2
    k2 = 2*np.pi*f/b2
    costheta2 = np.ndarray(len(theta2), dtype=complex)
    for u in range(len(theta2)):
        if m.isnan(theta2[u]):
            costheta2[u] = 1j*np.sqrt(((c2/c1*np.sin(theta[u]))**2-1))
        else:
            costheta2[u] = np.cos(theta2[u])

    cosgamma2 = np.ndarray(len(gamma2), dtype=complex)
    for u in range(len(gamma2)):
        if m.isnan(gamma2[u]):
            cosgamma2[u] = 1j*np.sqrt(((b2/c1*np.sin(theta[u]))**2-1))
        else:
            cosgamma2[u] = np.cos(gamma2[u])

    cos2gamma2 = np.ndarray(len(2*gamma2), dtype=complex)
    for u in range(len(2*gamma2)):
        if m.isnan(gamma2[u]):
            cos2gamma2[u] = 2*cosgamma2[u]**2-1
        else:
            cos2gamma2[u] = np.cos(2*gamma2[u])

    sin2gamma2 = np.ndarray(len(2*gamma2), dtype=complex)
    for u in range(len(2*gamma2)):
        if m.isnan(gamma2[u]):
            singamma2 = b2/c1*np.sin(theta[u])
            sin2gamma2[u] = 2*singamma2*cosgamma2[u]
        else:
            sin2gamma2[u] = np.sin(2*gamma2[u])

    P = K2*h*costheta2 # <--- x2
    Q = k2*h*cosgamma2 # <--- x2
    Z1 = rho1*c1/np.cos(theta1)
    Z2 = rho2*c2/costheta2
    Z2t = rho2*b2/cosgamma2
    M = (Z2/Z1)*cos2gamma2**2/np.tan(P) + (Z2t/Z1)*sin2gamma2**2/np.tan(Q)
    N = (Z2/Z1)*cos2gamma2**2/np.sin(P) + (Z2t/Z1)*sin2gamma2**2/np.sin(Q)
    A = 1j*(M**2-N**2-1)
    B = 1j*(M**2-N**2+1)

    T = 2*N/(2*M + A)
    R = B/(2*M + A)
    R_ana = R.copy()
    T_ana = T.copy()

    x=np.degrees(theta)
    theta_ana = x.copy()
    plt.plot(x, np.abs(R)**2, linestyle='--', label='$|R_{LL}|^2~(th)$')
    plt.plot(x, np.abs(T)**2, linestyle='--', label='$|T_{LL}|^2~(th)$')
    plt.plot(x, np.angle(R)/np.pi, linestyle='--')
    plt.plot(x, np.angle(T)/np.pi, linestyle='--')
    plt.legend()
    plt.show()

    # return (theta_num, R_num, T_num, theta_ana, R_ana, T_ana)
    return
