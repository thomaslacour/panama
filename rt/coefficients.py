import numpy as np
import math as m

import material.db as mat
import rt.waves as waves

"""
freq_vec = [freq1, freq2, ...]
layer1 = ('name', thickness, 'mat', 'inc', d, phi, poly)
stack = [ layer0, layer1, ... ]
stacks = [ stack(freq1), stack(freq2), ... ]
"""

def layer_properties(freq_vec, material):
    """
    Initialize and format layer properties depending on material
    type (meta or not meta) and return:

        d, cL, cT, rhoL, rhoT, alphaL, alphaT

    where each quantities are array.
    """
    # name of the material
    material_name = material[0]
    # thickness of the material (reshape with freq shape, in a tuple, to
    # allow the sum with the tuple of material properties)
    thickness = (np.array( [material[1]]*len(freq_vec) ), )
    # check if we have to pass extra arguments for non homogenous material
    if material_name == 'meta':
        param = material[2:]
    else:
        param = ()
    # read/compute material properties
    prop = mat.properties(material_name, freq_vec, *param)

    return thickness + prop




def norm_s(omega, c, alpha=0):
    """ modulus of the slowness vector """
    if omega !=0:
        return 1/c + 1j**alpha/omega
    else:
        return 0




def rt(freq_vec, angle_vec, stacks_list, pola='LL'):
    """ compute reflexion and transmission coefficients.  """

    cells = 1
    Coeff = 1

    R = np.ndarray( (len(freq_vec), len(angle_vec)), dtype=float )
    T = np.ndarray( (len(freq_vec), len(angle_vec)), dtype=float )

    for j, angle in enumerate(angle_vec):
        for i, (freq, stack) in enumerate(zip(freq_vec, stacks_list)):

            omega          = 2*np.pi*freq
            d              = stack[0][1:-1] # we exclude halfspaces
            cL, cT         = (*stack[1],), (*stack[2],)
            rhoL, rhoT     = (*stack[3],), (*stack[4],)
            alphaL, alphaT = (*stack[5],), (*stack[6],)

            lp = [d, cL, cT, rhoL, rhoT, alphaL, alphaT]

            # norm of the slowness vector
            if pola[0] == 'L':
                s = norm_s(omega, cL[0], alphaL[0])
            elif pola[0] == 'T':
                s = norm_s(omega, cT[0], alphaT[0])

            # slowness component
            s1 = s*m.sin(m.radians(angle))
            # scattering matrix
            S = waves.scattering_matrix(omega, s1, cells, *lp)

            if pola == 'LL':
                # displacement coefficients
                r = S[S.shape[0]//2,0] # rLL
                t = S[0,0]             # tLL
                # conversion to energy coefficients
                Coeff = rhoL[-1]*cL[-1]/(rhoL[0]*cL[0])

            R[i,j] = np.abs(r)**2       # RLL
            T[i,j] = Coeff*np.abs(t)**2 # TLL

    if len(freq_vec) == 1 or len(angle_vec) == 1:
        R = R.flat[:]
        T = T.flat[:]

    return R, T





def stacks(freq_vec, stack_config):
    ''' make a stack of layers properties '''
    # ----- stacks[layers][properties][frequency]
    stacks = [ layer_properties(freq_vec, l) for l in stack_config ]
    # ----- reorder the axis to get a list of the stacks:
    #       stacks[frequency][properties,layers]
    return list(np.array(stacks).transpose(2,1,0))





if __name__ == "__main__":

    lay = []

    # initialize the representation vectors
    freq_vec = np.linspace(0, 0.02, 101)
    theta_vec = [0]

    # configuration of the layers
    # None thickness for halfspaces
    left = ('water', None)
    lay += [ ('meta', 100, 'PU', 'billes', 0.25, 6.67, 10) ]
    # lay += [ (100, 'meta', 'PU', 'billes', 0.25, 6.67, 10) ]
    lay += [ ('steel', 100) ]
    right = ('air', None)
    stack_config = [left, *lay, right]

    stack_list = stacks(freq_vec, stack_config)
    R, T = rt(freq_vec, theta_vec, stack_list , pola='LL')

    import rt.display as plot_rt
    import matplotlib.pyplot as plt

    plot_rt.plot1d(freq_vec, R, T)
    # plot_rt.plot2d(freq_vec, theta_vec, R, T)

    # save layer properties
    lp = [ tuple(u.tolist()) for u in stack_list[0] ]
    lp = lp[0][1:-1], lp[1:]

    plt.show(block=False)
