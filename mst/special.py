import scipy.special as sp
from numpy import array, ndarray

# -- polynome de legendre -- https://stackoverflow.com/a/45775310
def legendre(n,X) :
    res = []
    for m in range(n+1):
        res.append(sp.lpmv(m,n,X))
    return res

# -- fonction de hankel sphérique --
# du premier premier type d'ordre de n pour une valeur de x
# notation : h1(n,x)
def h1(n,z):
    return sp.spherical_jn(n,z) + 1j*sp.spherical_yn(n,z)
# sa dérivé
def h1p(n,z):
    return n*h1(n,z)/z - h1(n+1,z)

# -- fonction de bessel sphérique --
# d'ordre de n pour une valeur de x
# notation : j(n,x)
def besseljspher(n,z):
    return sp.spherical_jn(n,z)
def jn(n,z):
    return sp.spherical_jn(n,z)
def jnp(n,z):
    # return n*jn(n,z)/z-jn(n+1,z)
    return sp.spherical_jn(n, z, derivative=True)
