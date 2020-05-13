import numpy as np
import numpy.linalg as npl

# direction 1: direction parallel to layers (Snell laws)
# direction 2: direction orthogonal to layers

#TODO: attention à externaliser et à sortir les propriétés des milieux
#      exterieurs, sinon fonction pas bien REUTILISABLE en dehors 
#      cas présent ... !

def s3_of_s1(s1, c, a=0):
    """
    SLOWNESS COMPONENT NORMAL TO THE SURFACE
    s1 (k1/omega): slowness component parallel to layers (Snell component)
    c:  phase speed of the wave
    a: (alpha/omega): attenuation of the wave
    """
    s = 1/c + 1j*a
    s3_2 = s**2 - s1**2
    s3 = np.sqrt(s3_2)
    if np.imag(s3)<0:
        s3 = -s3
    return s3


def wave_polar(s1, s3, polar):
    """
    POLARISATION VECTOR
    s1 (k1/omega) and s3 (k3/omega): slowness components of partial wave
    polar: 'L' for longitudinal wave, 'T' for transverse wave
    """
    norm = np.sqrt(s1**2+s3**2)
    if polar == 'L':
        P = np.array([s1, s3])/norm
    elif polar == 'T':
        P = np.array([-s3, s1])/norm
    return P


def transfert_matrix_1layer(omega, s1, d, CL, CT, rhoL, rhoT, alphaL, alphaT):
    """
    TRANSFERT MATRIX FOR 1 LAYER
    d, CL, CT, rhoL, rhoT are scalars defining the material properties of the layer
    Xi: matrix 4-coulumn of state-vectors of each partial wave
    phase_matrix: 4x4 diagonal matrix of phase factors
    """
    list_pola, ind_st_vec = partial_wave('solid')
    c_     = {'L':CL, 'T':CT}
    rho_   = {'L':rhoL, 'T':rhoT}
    alpha_ = {'L':alphaL/omega, 'T':alphaT/omega}
    s3     = np.array([ s3_of_s1(s1, c_[i], alpha_[i]) for i in list_pola ])
    s3vec  = np.append( s3, -s3 )
    s_     = [s1, s3]
    Xi = Xi_matrix(omega, s_, c_, rho_, list_pola)[ind_st_vec,:]
    phase_matrix = np.diag(np.exp(1j*omega*s3vec*d))

    return Xi @ phase_matrix @ npl.inv(Xi)


#TODO: on vire les milieux externes le test doit être fait en dehors !
def transfert_matrix(omega, s1, unit_cells, d, CL, CT, rhoL, rhoT, alphaL, alphaT):
    """
    TRANSFERT MATRIX
    d, CL, CT, rhoL, rhoT, alphaL, alphaT are arrays of length
    the number of layers. It return the transfert matrix M:
       (4x4) if outer media are SOLIDS
       (2x2) if outer media are FLUIDS
    """
    M = np.eye(4,dtype=complex)
    for i in range(len(d)):
        # M = transfert_matrix_1layer(omega, s1, d[i], CL[1:-1][i], CT[1:-1][i], rhoL[1:-1][i], rhoT[1:-1][i]) @ M
        C = (CL[1:-1][i], CT[1:-1][i])
        rho = (rhoL[1:-1][i], rhoT[1:-1][i])
        alpha = (alphaL[1:-1][i], alphaT[1:-1][i])
        M = transfert_matrix_1layer(omega, s1, d[i], *C, *rho, *alpha) @ M
    if unit_cells > 1:
        M = npl.matrix_power(M, unit_cells)
    if outer_material(test=CT)['fluid']:
        m11 = -npl.det( M[1:3,0:2] )/M[2,0]
        m12 = npl.det( M[1:3,(-1,0)] )/M[2,0]
        m21 = npl.det( M[2:4,0:2] )/M[2,0]
        m22 = -npl.det( M[2:4,(-1,0)] )/M[2,0]
        M = np.array([[m11,m12],
                      [m21,m22]])
    return M

# le test du milieu est à l'extérieur de la fonction
# def reduced_transfert_matrix(omega, s1, unit_cells, d, CL, CT, rhoL, rhoT, alphaL, alphaT)
#     M = transfert_matrix(omega, s1, unit_cells, d, CL, CT, rhoL, rhoT, alphaL, alphaT)
#     m11 = -npl.det( M[1:3,0:2] )/M[2,0]
#     m12 = npl.det( M[1:3,(-1,0)] )/M[2,0]
#     m21 = npl.det( M[2:4,0:2] )/M[2,0]
#     m22 = -npl.det( M[2:4,(-1,0)] )/M[2,0]
#     M = np.array([[m11,m12],
#                   [m21,m22]])
#     return M

# def Xi_vec_solid(omega, s, c_, rho_, polar):
def Xi_vec(omega, s, c_, rho_, polar):
    """
    STATE VECTOR (1 partial wave)
    s is a tuple with the slowness component of the partial wave (s1, s3)
    c_, rho_ are dicts with 2 arrays {L:, T:} of length the number of layers
    polar is a list of polarisation of propagating waves
    """
    P = wave_polar(*s, polar)
    s1, s3 = s[0], s[1]
    p1, p3 = P[0], P[1]
    v1, v3 = -p1, -p3
    stress13 = rho_['T']*(c_['T']**2)*( s1*p3 + s3*p1 )
    stress33 = rho_['L']*(c_['L']**2)*( s3*p3 ) \
              + (rho_['L']*c_['L']**2-2*rho_['T']*c_['T']**2)*( s1*p1 )
    eta = 1j*omega*np.array([v1, v3, stress13, stress33])
    return eta

def Xi_vec_fluid(omega, s, c_, rho_):
    """
    STATE VECTOR (1 partial wave)
    s is a tuple with the slowness component of the partial wave (s1, s3)
    c_, rho_ are dicts with 2 arrays {L:, T:} of length the number of layers
    polar is a list of polarisation of propagating waves
    """
    P = wave_polar(*s, polar='L')
    s1, s3 = s[0], s[1]
    p1, p3 = P[0], P[1]
    v1, v3 = -p1, -p3
    minus_p = rho_['L']*(c_['L']**2)*( s3*p3 + s1*p1)
    eta = 1j*omega*np.array([ v3, minus_p ])
    return eta


# def Xi_matrix_solid(omega, s_, c_, rho_):
def Xi_matrix(omega, s_, c_, rho_, list_pola):
    """
    STATE VECTOR (all partial waves)
    s_ is a list with [ s1, [s3L, s3T, -s3L, -s3T] ]
    c_, rho_ are dicts with 2 arrays {L:, T:} of length the number of layers
    list_pola is a list of all propagating waves
    """
    # pola_ = ['L','T']*2
    pola_ = list_pola*2
    N = int(len(pola_))
    s1, s3 = s_[0], s_[1]
    # s3    = np.array([ s3_of_s1(s1, c_[i], alpha_[i]) for i in list_pola ])
    s3vec = np.append( s3, -s3 )
    Xi = np.array([ Xi_vec(omega, (s1, s3vec[j]), c_, rho_, pola_[j]) for j in range(N)]).T
    return Xi


def outer_material(**kwargs):
    """
    TEST FOR FLUID OR SOLID
    """
    CT = kwargs.get('test', None)
    if CT[0]==0 or CT[-1]==0:
        fluid, solid = True, True
    if CT[0]==0 and CT[-1]==0:
        fluid, solid = True, False
    if CT[0]!=0 and CT[-1]!=0:
        solid, fluid = True, False

    return { 'fluid':fluid, 'solid':solid }


def partial_wave(type):
    """
    LIST OF PROPAGATING WAVES
    type: solid or fluid
    """
    if type == 'solid':
        list = ['L', 'T']
        ind  = (0,1,2,3)
    if type == 'fluid':
        list = ['L']
        ind  = (1,3)

    return list, ind


# def scattering_matrix()


def scattering_matrix(omega, s1, unit_cells, d, CL, CT, rhoL, rhoT=None, alphaL=None, alphaT=None):
    """
    SCATTERING MATRIX relating converging waves to diverging ones.
    s1 (k1/omega): slowness component parallel to layers (Snell component)
    d, CL, CT, rhoL, rhoT are scalars defining the material properties of the layer
    """

    if rhoT == None: rhoT = rhoL
    if alphaL == None: alphaL = [0]*len(CL)
    if alphaT == None: alphaT = [0]*len(CT)

    for i in ['fluid', 'solid']:
        if outer_material(test=CT)[i]:
            list_pola, ind_st_vec = partial_wave(i)

    # CL, CT sans prop ext
    # if gauche=solid et droit=solid
    #   transfertmatrix_44
    # elif gauche=fluid et droit=fluid
    #   transfertmatrix_22
    # elif f/s or s/f
    #   ...
    M = transfert_matrix(omega, s1, unit_cells, d, CL, CT, rhoL, rhoT, alphaL, alphaT)

    N = int(len(list_pola))
    pola_ = list_pola*2

    Xi = []

    # TODO: sortir le calcul du vecteur d'état ?
    for i in (0, -1):
        c_     = {'L':CL[i], 'T':CT[i]}
        alpha_ = {'L':CL[i], 'T':CT[i]}
        rho_   = {'L':rhoL[i], 'T':rhoT[i]}
        alpha_ = {'L':alphaL[i]/omega, 'T':alphaT[i]/omega}
        s3     = np.array([ s3_of_s1(s1, c_[i], alpha_[i]) for i in list_pola ])
        s3_vec = np.append( s3, -s3 )
        s_     = [s1, s3]
        Xi.append(
            Xi_matrix(omega, s_, c_, rho_, list_pola )[ind_st_vec,:]
        )

    if npl.det(M @ Xi[0])==0:
        B = np.eye(N,N)
    else:
        B = npl.inv( M @ Xi[0] ) @ Xi[1]
    B1 = B[:N,:N]
    B2 = B[:N,N:]
    B3 = B[N:,:N]
    B4 = B[N:,N:]
    S1 = npl.inv(B1)
    S2 = -S1 @ B2
    S3 = B3 @ S1
    S4 = B4 - S3 @ B2

    return np.block([[S1,S2],[S3,S4]])
