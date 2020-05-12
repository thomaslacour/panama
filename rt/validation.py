import numpy as np
import cmath as m
import matplotlib.pyplot as plt

from rt import actions as action
from rt import display as display

from material import db as material

def test_1layer_NormalIncidence(*arg):
    print('\n1| Normal incidence for 1 layer immersed in a fluid')
    # --------------------------------------------------------
    # Test de validation du code en incidence normale.
    # --------------------------------------------------------
    fmin, fmax, fnum = 0.01, 10, 301
    data = arg[0]
    data['unit_cells']      = 1
    data['layers']['1']     = ('steel', 1.0)
    data['halfspace_left']  = 'water'
    data['halfspace_right'] = 'water'
    data['f_min']           = fmin
    data['f_max']           = fmax
    data['f_num']           = fnum
    data['theta_fix']       = 0 # incidence normale
    data['todo']            = 0 # RT_vs_frequency
    # layers_prop = action.extract_layers_properties(data)
    lp = display.init_lp(data, omega=0)
    x, R, T, ax1 = display.rt_panama(data, plot_vs='frequency', phase=True)
    f_num, R_num, T_num = x, R, T #<--------
    # --------------------------------------------------------
    # formule analytique
    # --------------------------------------------------------
    # ARISTÉGUI, et al., Wave Motion (2010), Eq. (2)
    # --------------------------------------------------------
    h = data['layers']['1'][1]/2
    freqs = np.linspace(fmin,fmax,fnum)
    rho1 = lp[3][1]
    cL1 = lp[1][1]
    Z1 = rho1*cL1
    rho0 = lp[3][0]
    cL0 = lp[1][0]
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

    ax1[0].plot(freqs, np.abs(R)**2, linestyle='--', label='$|R_{LL}|^2~(th)$')
    ax1[0].plot(freqs, np.abs(T)**2, linestyle='--', label='$|T_{LL}|^2~(th)$')
    ax1[1].plot(freqs, np.angle(R)/np.pi, linestyle='--')
    ax1[1].plot(freqs, np.angle(T)/np.pi, linestyle='--')
    ax1[0].legend()
    # kd
    ax2 = ax1[0].twiny()
    f = np.array([0,0.5,1])*cL1/(2*h)
    def tick_function(freq):
        x = freq*2*h/cL1
        return ["%.3f" % z for z in x]
    ax2.set_xlim(ax1[0].get_xlim())
    ax2.set_xticks(f)
    ax2.set_xticklabels(tick_function(f))
    ax2.set_xlabel("$d\,/\,\lambda$")
    plt.show()

    return (f_num*2*h/cL1, R_num, T_num, f_ana*2*h/cL1, R_ana, T_ana)

def test_1layer_Angle(*arg):
    print('\n2| Incidence for 1 layer in a fluid')
    # --------------------------------------------------------
    # Test de validation du code sous incidence.
    # --------------------------------------------------------
    thetamin, thetamax, thetanum = 0, 90, 901
    f = 1 # MHz
    data                    = arg[0]
    data['unit_cells']      = 1
    data['layers']['1']     = ('steel', 2.0)
    data['halfspace_left']  = 'water'
    data['halfspace_right'] = 'water'
    data['theta_min']       = thetamin
    data['theta_max']       = thetamax
    data['theta_num']       = thetanum
    data['f_fix']           = f
    data['todo']            = 1 # RT_vs_angle
    # layers_prop = action.extract_layers_properties(data, lm)
    lp = display.init_lp(data, omega=0)
    x, R, T, ax = display.rt_panama(*arg, plot_vs='angle', phase=True)
    theta_num = x.copy()
    R_num = R.copy()
    T_num = T.copy()

    # formule analytique: 3(fluid) | 2(solid) | 1(fluid=3)
    # --------------------------------------------------------
    # BREKHOVSKIKH, Waves in Layered Media, p. 71
    # Eq. (10.7) et (10.8)
    # --------------------------------------------------------
    h = data['layers']['1'][1]
    theta = np.radians(np.linspace(thetamin,thetamax,thetanum))
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
    ax[0].plot(x, np.abs(R)**2, linestyle='--', label='$|R_{LL}|^2~(th)$')
    ax[0].plot(x, np.abs(T)**2, linestyle='--', label='$|T_{LL}|^2~(th)$')
    ax[1].plot(x, np.angle(R)/np.pi, linestyle='--')
    ax[1].plot(x, np.angle(T)/np.pi, linestyle='--')
    ax[0].legend()
    plt.show()

    return (theta_num, R_num, T_num, theta_ana, R_ana, T_ana)
