from numpy import ndarray, array, arange, zeros, ones
from numpy import transpose, concatenate, hstack, vstack,  reshape
from numpy import cos, sin, sqrt, pi
from numpy.linalg import det, solve
from .special import *
from .special import h1, besseljspher
from numpy import sum as array_sum

def SnLL(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax):
    if all(mu0)==0 and all(mu1)==0: # ---------------- matrice fluide / inclusions fluides
        # print('fluid/fluid')
        fctCalculScatteringCoeffs = ffLL

    elif all(mu0)==0 and any(mu1)!=0: # --------------- matrice fluide / inclusions solides
        # print('fluid/solid')
        fctCalculScatteringCoeffs = fsLL

    elif any(mu0)!=0 and all(mu1)==0: # ----------------- matrice solide / inclusions fluides
        # print('solid/fluid')
        fctCalculScatteringCoeffs = sfLL

    elif any(mu0)!=0 and any(mu1)!=0: # ----------------- matrice solide / inclusions solides
        # print('solid/solid')
        fctCalculScatteringCoeffs = ssLL_new

    S  = fctCalculScatteringCoeffs(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)

    Sn = S[-1,0]

    # test convergence (à reprendre)
    # [Sn,S, nmax]= testconvergence(fctCalculScatteringCoeffs, lambda0,mu0,rho0,lambda2,mu2,rho2,f,rd,nmax);

    # ------ subfunctions
    # -------------------------------------------
    # function [Sn, S, nmax]= testconvergence(fctCalculScatteringCoeffs, lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)
    # % fonction de jerome
    # %test de convergence de la fonction de forme (vérifie que le dernier coeff
    # %est 10^(-10) fois plus petit que le 1er)
    # 
    # % nmax nbre de coeffs; S(1:n+1) = S0, S1,..Sn
    # temp=1;
    #  while temp>10^(-10)
    #             S=fctCalculScatteringCoeffs(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax);
    #             Sn=S(end,1);
    #             for cons=2:nmax+1
    #                 Sn=Sn+2*S(end,cons);
    #             end
    #             temp=abs(S(end,nmax+1))/abs(Sn);
    #             if temp>10^(-10)
    #                 nmax=nmax+10;
    #             end
    # end%test de convergence de la fonction de forme

    return S, nmax



# function [S] = ff(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,n_modes)
# %programme calculant les coefficients de diffusion modaux d'un système
# %fluide fluide en fonction de f vecteur (colonne) de la fréquence
# %notation : ff(lambda0,rho0,lambda1,rho1,f,rd,n_modes)
# %lambda et rho peuvent être des vecteurs colonne de même taille que f,
# %rd est le rayon du la sphère et n_modes est le mode maximal que l'on 
# %souhaite calculer
def ffLL(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, n_modes):
    c0 = sqrt(lambda0/rho0)
    c1 = sqrt(lambda1/rho1)
    k0a = 2*pi*f*rd/c0
    k1a = 2*pi*f*rd/c1
    X = k0a
    Y = k1a
    nfreq = len(f)

    S = zeros( (2,nfreq,n_modes+1), dtype=complex )

    for n in range(0,n_modes+1):
        A = zeros( (2,2,nfreq), dtype=complex )
        B = zeros( (2,nfreq), dtype=complex )
        A[0,0,:] = n*h1(n,X)  -  X*h1(n+1,X); #in C
        A[1,0,:] = -(X**2)*h1(n,X) #in C
        A[0,1,:] = Y*besseljspher(n+1,Y) - n*besseljspher(n,Y) #in R
        A[1,1,:] = lambda1/lambda0*(Y**2)*besseljspher(n,Y) #in R

        B[0]  = X*besseljspher(n+1,X) - n*besseljspher(n,X) #in R
        B[1]  = (X**2)*besseljspher(n,X) #in R
        for q in range(0,nfreq):
            S[:,q,n] = solve(A[:,:,q], B[:,q])

    if f[0]==0:
        S[0,:]=0;

    return S


# function [S] = sf(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)
# %programme calculant les coefficients de diffusion modaux d'un système
# %fluide fluide en fonction de f vecteur (colonne) de la fréquence
# %notation : sf(lambda1,mu1,rho1,lambda2,rho2,f,rd,nmax)
# %lambda et rho peuvent être des vecteurs colonne de même taille que f,
# %rd est le rayon du la sphère et nmax est le mode maximal que l'on 
# %souhaite calculer
def sfLL(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, nmax):
    c0 = sqrt((lambda0+2*mu0)/rho0)
    a = 2*pi*f*rd/c0
    b = a*((lambda0+2*mu0)/mu0)**(0.5)
    c = a*((lambda0+2*mu0)*rho1/(rho0*lambda1))**(0.5)
    nb = len(f)

    S = zeros( (3,nb,nmax+1), dtype=complex )

    for p in range(1,nmax+2):
        A = zeros( (3,3,nb), dtype=complex )
        B = zeros( (3,nb), dtype=complex )
        n = p-1
        A[0,0,:] = a*h1(n+1,a)-n*h1(n,a)
        A[1,0,:] = 4*a*h1(n+1,a)-(b**2+2*n*(1-n))*h1(n,a)
        A[2,0,:] = 2*(1-n)*h1(n,a)+2*a*h1(n+1,a)
        A[0,1,:] = n*(1+n)*h1(n,b)
        A[1,1,:] = 2*n*(1+n)*((1-n)*h1(n,b)+b*h1(n+1,b))
        A[2,1,:] = (2*(n**2-1)-b**2)*h1(n,b)+2*b*h1(n+1,b)
        A[0,2,:] = n*besseljspher(n,c)-c*besseljspher(n+1,c)
        A[1,2,:] = lambda1/mu0*c**2*besseljspher(n,c)
        B[0,:] = -a*besseljspher(n+1,a)+n*besseljspher(n,a)
        B[1,:] = (b**2+2*n*(1-n))*besseljspher(n,a)-4*a*besseljspher(n+1,a)
        B[2,:] = 2*(n-1)*besseljspher(n,a)-2*a*besseljspher(n+1,a)
        for q in range(0,nb):
            S[0:,q,p-1] = solve(A[:,:,q], B[:,q])
            #[TnLL, TnTL, CnLL]

    if f[0]==0: 
        S[:,0,:] = 0

    return S


# TODO: A REECRIRE PLUS PROPREMENT
# function [S] = ss(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)
# %programme calculant les coefficients de diffusion modaux d'un système
# %fluide fluide en fonction de f vecteur (colonne) de la fréquence
# %notation : solsol(lambda1,mu1,rho1,lambda2,mu2,rho2,f,rd,nmax)
# %lambda et rho peuvent être des vecteurs colonne de même taille que f,
# %rd est le rayon du la sphère et nmax est le mode maximal que l'on 
# %souhaite calculer
def ssLL(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, nmax):
    c0 = sqrt((lambda0+2*mu0)/rho0)
    a = 2*pi*f*rd/c0
    b = a*((lambda0+2*mu0)/mu0)**(0.5)
    c = a*((lambda0+2*mu0)*rho1/(rho0*(lambda1+2*mu1)))**(0.5)
    d = a*((lambda0+2*mu0)*rho1/(rho0*mu1))**(0.5)
    nb = len(f)
    S = zeros( (nb,nmax+1), dtype=complex )

    for p in range(1,nmax+2):
        n = p-1
        a11 = n*h1(n,a)-a*h1(n+1,a)
        a21 = h1(n,a)
        a31 = (2*n*(n-1)-b**2)*h1(n,a)+4*a*h1(n+1,a)
        a41 = 2*((n-1)*h1(n,a)-a*h1(n+1,a))
        a12 = -n*(n+1)*h1(n,b)
        a22 = b*h1(n+1,b)-(n+1)*h1(n,b)
        a32 = 2*n*(n+1)*((1-n)*h1(n,b)+b*h1(n+1,b))
        a42 = (2*(1-n**2)+b**2)*h1(n,b)-2*b*h1(n+1,b)
        a13 = c*besseljspher(n+1,c)-n*besseljspher(n,c)
        a23 = -besseljspher(n,c)
        a33 = -mu1/mu0*((2*n*(n-1)-d**2)*besseljspher(n,c)+4*c*besseljspher(n+1,c))
        a43 = -2*mu1/mu0*((n-1)*besseljspher(n,c)-c*besseljspher(n+1,c))
        a14 = n*(n+1)*besseljspher(n,d)
        a24 = (n+1)*besseljspher(n,d)-d*besseljspher(n+1,d)
        a34 = -mu1/mu0*2*n*(n+1)*((1-n)*besseljspher(n,d)+d*besseljspher(n+1,d))
        a44 = -mu1/mu0*((2*(1-n^2)+d**2)*besseljspher(n,d)-2*d*besseljspher(n+1,d))
        b1 = a*besseljspher(n+1,a)-n*besseljspher(n,a);
        b2 = -besseljspher(n,a)
        b3 = -(2*n*(n-1)-b**2)*besseljspher(n,a)-4*a*besseljspher(n+1,a)
        b4 = -2*((n-1)*besseljspher(n,a)-a*besseljspher(n+1,a))
        for q in range(0,nb):
            A = zeros( (4,4), dtype=complex )
            A[0,0] = a11[q]
            A[0,1] = a12[q]
            A[0,2] = a13[q]
            A[0,3] = a14[q]
            A[1,0] = a21[q]
            A[1,1] = a22[q]
            A[1,2] = a23[q]
            A[1,3] = a24[q]
            A[2,0] = a31[q]
            A[2,1] = a32[q]
            A[2,2] = a33[q]
            A[2,3] = a34[q]
            A[3,0] = a41[q]
            A[3,1] = a42[q]
            A[3,2] = a43[q]
            A[3,3] = a44[q]
            B = zeros( (4,4), dtype=complex )
            B[0,0] = b1[q]
            B[0,1] = a12[q]
            B[0,2] = a13[q]
            B[0,3] = a14[q]
            B[1,0] = b2[q]
            B[1,1] = a22[q]
            B[1,2] = a23[q]
            B[1,3] = a24[q]
            B[2,0] = b3[q]
            B[2,1] = a32[q]
            B[2,2] = a33[q]
            B[2,3] = a34[q]
            B[3,0] = b4[q]
            B[3,1] = a42[q]
            B[3,2] = a43[q]
            B[3,3] = a44[q]

            S[q,p-1] =det(B)/det(A)

    if f[0]==0:
        S[1,:] = 0

    return S


# function [S] = fs(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)
# %programme calculant les coefficients de diffusion modaux d'un système
# %fluide fluide en fonction de f vecteur (colonne) de la fréquence
# %notation : fs(lambda0,rho0,lambda1,mu1,rho1,f,rd,nmax)
# %lambda et rho peuvent être des vecteurs colonne de même taille que f,
# %rd est le rayon du la sphère et nmax est le mode maximal que l'on 
# %souhaite calculer
def fsLL(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, nmax):
    c0 = sqrt(lambda0/rho0)
    a = 2*pi*f*rd/c0
    b = a*((lambda0*rho1)/(rho0*(lambda1+2*mu1)))**(0.5);
    c = a*(lambda0*rho1/(rho0*mu1))**(0.5)
    nFreq = len(f)
    S = zeros( (3,nFreq,nmax+1), dtype=complex )

    for p in range(1,nmax+2):
        A = zeros( (3,3,nFreq), dtype=complex )
        B = zeros( (3,nFreq), dtype=complex )
        n = p-1

        A[0,0] = n*h1(n,a)-a*h1(n+1,a)
        A[1,0] = a**2*h1(n,a)
        A[0,1] = b*besseljspher(n+1,b)-n*besseljspher(n,b)
        A[1,1] = mu1/lambda0*(4*b*besseljspher(n+1,b)-(c**2+2*n*(1-n))*besseljspher(n,b))
        A[2,1] = 2*(1-n)*besseljspher(n,b)+2*b*besseljspher(n+1,b)
        A[0,2] = n*(1+n)*besseljspher(n,c)
        A[1,2] = mu1/lambda0*2*n*(1+n)*((1-n)*besseljspher(n,c)+c*besseljspher(n+1,c))
        A[2,2] = (2*(n^2-1)-c**2)*besseljspher(n,c)+2*c*besseljspher(n+1,c)

        B[0]=a*besseljspher(n+1,a)-n*besseljspher(n,a)
        B[1]=-a**2*besseljspher(n,a)
        for q in range(0,nFreq):
            S[0:,q,p-1] = solve(A[:,:,q], B[:,q])
            #[TnLL, TnTL, CnLL]

    if f[0]==0:
        S[1,:]=0

    return S

def calcul_f_sc(theta,f,k0, S, nCoeffs):
    n = arange(0,nCoeffs+1,1)
    TnLL = S[:,0:nCoeffs+1]
    if theta==0:
        return array_sum( (2*n+1)*TnLL, axis=1)/(1j*k0)
    elif theta==pi:
        return array_sum( (2*n+1)*((-1)**n)*TnLL, axis=1)/(1j*k0)
    nfreq = len(f)
    Pncos0 = calculPn(cos(theta),nCoeffs);
    Pncos0 = ones( (nfreq,1) )*transpose(Pncos0)
    n = ones( (nfreq,1) )*n;
    fsc_theta_f = array_sum( (2*n+1)*TnLL*Pncos0, axis=1)
    fsc_theta_f = 1/(1j*k0)*fsc_theta_f

    return fsc_theta_f


def calculPn(x,nCoeffs):
    # P0, P1,...Pnmax-1
    Pn = zeros( (nCoeffs +1,1) )

    for n in range(0,nCoeffs+1):
        p = legendre(n, x );
        Pn[n] = p[0]
    return Pn


##-----------
def ssLL_new(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, nmax):
    c0L = ((lambda0+2*mu0)/rho0)**0.5
    c0T = (mu0/rho0)**0.5
    c1L = ((lambda1+2*mu1)/rho1)**0.5
    c1T = (mu1/rho1)**0.5
    a, b, c, d = [ 2*pi*f*rd/i for i in [c0L, c0T, c1L, c1T ] ]
    nb = len(f)
    S = zeros( (4,nb,nmax+1), dtype=complex )

    for p in range(1,nmax+2):
        A = zeros( (4,4,nb), dtype=complex )
        B = zeros( (4,nb), dtype=complex )
        n = p-1
        A[0,0,:] = a*h1p(n,a)
        A[0,1,:] = n*(n+1)*h1(n,b)
        A[0,2,:] = -c*jnp(n,c)
        A[0,3,:] = -n*(n+1)*jn(n,d)
        A[1,0,:] = h1(n,a)
        A[1,1,:] = (b*h1p(n,b) + h1(n,b))
        A[1,2,:] = -jn(n,c)
        A[1,3,:] = -d*jnp(n,d) - jn(n,d)
        #          o---> !! diff Lepert(-)/Ba(+)
        #          |
        A[2,0,:] = -((b**2 - 2*n*(n+1))*h1(n,a) + 4*a*h1p(n,a))
        A[2,1,:] = 2*n*(n+1)*(b*h1p(n,b) - h1(n,b))
        A[2,2,:] = mu1/mu0*((d**2-2*n*(n+1))*jn(n,c) + 4*c*jnp(n,c))
        A[2,3,:] = 2*n*(n+1)*mu1/mu0*(-d*jnp(n,d)+jn(n,d))
        A[3,0,:] = (2*(h1(n,a)-a*h1p(n,a)))
        A[3,1,:] = (b**2 - 2*(n**2+n-1))*h1(n,b) + 2*b*h1p(n,b)
        A[3,2,:] = 2*mu1/mu0*(c*jnp(n,c) - jn(n,c))
        A[3,3,:] = -mu1/mu0*((d**2-2*(n**2+n-1))*jn(n,d) + 2*d*jnp(n,d))
        B[0,:] = -a*jnp(n,a)
        B[1,:] = -jn(n,a)
        #                o---> !! diff Lepert()/Ba(n)
        #                |
        B[2,:] = (b**2-2*n*(n+1))*jn(n,a) + 4*a*jnp(n,a)
        B[3,:] = 2*(a*jnp(n,a) - jn(n,a))
        for q in range(0,nb):
            S[:,q,p-1] = solve(A[:,:,q], B[:,q])
            #[TnLL, TnTL, CnLL, DnTL]

    if f[0]==0: S[1,:] = 0

    return S

def SnTT(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax):
    if all(mu0)==0 and all(mu1)==0:
        fctCalculScatteringCoeffs = ffTT
    elif all(mu0)==0 and any(mu1)!=0:
        fctCalculScatteringCoeffs = fsTT
    elif any(mu0)!=0 and all(mu1)==0:
        fctCalculScatteringCoeffs = ssTT
    elif any(mu0)!=0 and any(mu1)!=0:
        fctCalculScatteringCoeffs = ssTT

    S  = fctCalculScatteringCoeffs(lambda0,mu0,rho0,lambda1,mu1,rho1,f,rd,nmax)
    return S, nmax

def ssTT(lambda0, mu0, rho0, lambda1, mu1, rho1, f, rd, nmax):
    c0L = ((lambda0+2*mu0)/rho0)**0.5
    c0T = (mu0/rho0)**0.5
    c1L = ((lambda1+2*mu1)/rho1)**0.5
    c1T = (mu1/rho1)**0.5
    if (c0T==0).any():
        c0T[:]=0.000001
        # print('warning: velocity of the matrice is zero, which may lead to divide by zero.\n' + '         We force it''s value to {}' . format(c0T[0]) )
    if (c1T==0).any():
        c1T[:]=0.000001
        # print('warning: velocity of the inclusion is zero, which may lead to divide by zero.\n' + '         We force it''s value to {}' . format(c1T[0]) )
    a, b, c, d = [ 2*pi*f*rd/i for i in [c0L, c0T, c1L, c1T ] ]
    nb = len(f)
    S = zeros( (6,nb,nmax+1), dtype=complex )

    for p in range(1,nmax+2):
        A = zeros( (4,4,nb), dtype=complex )
        B = zeros( (4,nb), dtype=complex )
        n = p-1
        A[0,0,:] = a*h1p(n,a)
        A[0,1,:] = n*(n+1)*h1(n,b)
        A[0,2,:] = -c*jnp(n,c)
        A[0,3,:] = -n*(n+1)*jn(n,d)
        A[1,0,:] = h1(n,a)
        A[1,1,:] = b*h1p(n,b) + h1(n,b)
        A[1,2,:] = -jn(n,c)
        A[1,3,:] = -(d*jnp(n,d) + jn(n,d))
        A[2,0,:] = ((b**2 - 2*n*(n+1))*h1(n,a) + 4*a*h1p(n,a))
        A[2,1,:] = 2*n*(n+1)*(h1(n,b)-b*h1p(n,b))
        A[2,2,:] = -mu1/mu0*( (d**2-2*n*(n+1))*jn(n,c) + 4*c*jnp(n,c) )
        A[2,3,:] = 2*n*(n+1)*mu1/mu0*( d*jnp(n,d) - jn(n,d) )
        A[3,0,:] = 2*(h1(n,a) - a*h1p(n,a))
        A[3,1,:] = ( (b**2 - 2*(n**2+n-1))*h1(n,b) + 2*b*h1p(n,b) )
        A[3,2,:] = 2*mu1/mu0*( c*jnp(n,c) - jn(n,c) )
        A[3,3,:] = -mu1/mu0*( (d**2-2*(n**2+n-1))*jn(n,d) + 2*d*jnp(n,d))
        B[0,:] = -n*(n+1)*jn(n,b)
        B[1,:] = -(b*jnp(n,b) + jn(n,b))
        B[2,:] = -(2*n*(n+1)*(jn(n,b)-b*jnp(n,b)))
        B[3,:] = -( (b**2 - 2*(n**2+n-1))*jn(n,b) + 2*b*jnp(n,b))
        for q in range(0,nb):
            S[0:4,q,p-1] = solve(A[:,:,q], B[:,q])
            #[TnTT, TnLT, CnLT, DnTT] <-- in plane (A. Ba)

    for p in range(1,nmax+2):
        A = zeros( (2,2,nb), dtype=complex )
        B = zeros( (2,nb), dtype=complex )
        n = p-1
        A[0,0,:] = -b*h1(n,b)
        A[0,1,:] = d*jn(n,d)
        A[1,0,:] = -b*( b*h1p(n,b) - h1(n,b) )
        A[1,1,:] = d*mu1/mu0*( d*jnp(n,d) - jn(n,d) )
        B[0,:] =  b*jn(n,b)
        B[1,:] = b*( b*jnp(n,b) - jn(n,b) )
        for q in range(0,nb):
            S[4:,q,p-1] = solve(A[:,:,q], B[:,q])
            #[tnTT, dnTT] <-- out of plane (A. Ba)

    if f[0]==0:
        # S(1,:)=0;
        S[1,:]=0
    # end

    return S

def f_tt(theta, f, k0T, S, nCoeffs):
    TnTT=S[0][:,1:nCoeffs+1]
    tnTT=S[1][:,1:nCoeffs+1]
    nfreq = len(f)
    n = arange(1, nCoeffs+1, 1)
    n = ones( (nfreq,1) )*n;
    if theta == 0:
        return array_sum( (2*n+1)*(TnTT+tnTT), axis=1 )/(1j*k0T)
    elif theta == pi:
        return array_sum( (2*n+1)*(TnTT-tnTT), axis=1 )/(1j*k0T)
