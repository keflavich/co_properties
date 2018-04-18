# CO LVG code
#Adam Ginsburg
# essentially solves eqn 16.3 of Rohlfs:
# n_l (C_lu + B_lu U) = n_u (A_ul + B_ul U + C_ul)

import sys
sys.output_line_width=180
from pylab import *
import scipy
import numpy
from numpy import append
#import LinearAlgebra
#from LinearAlgebra import solve_linear_equations
from scipy import linalg
from numpy import float64

#constants
h = 6.626068e-27    # Planck's Constant (erg s)
c = 2.99792458e10   # Speed of Light (cm/s)
k = 1.3806503e-16   # Boltzmann's Constant (erg / K)
m_co = 4.65e-23     # Mass of CO
m_h2 = 3.32e-24     # Mass of H_2
m_red = m_co*m_h2 / (m_co+m_h2)  # Reduced mass of CO/H_2 interaction
#cross section
sigma = 2e-15 #cm^-2

nlevels = 8 # how many levels to use?

#transitions: 1-0  2-1  3-2  4-3  5-4  6-5  7-6
levels = arange(nlevels)  # how many levels...
gi = levels*2+1.    # degeneracy
# the frequencies of the transitions are not given by the fibonacci sequences * 115.3e9 Hz
# fib = asarray([1,3,6,10,15,21,27,34],dtype='float64') * 115.3e9
Be = 57.63596828e9 # rotational constant for 12CO
#nuij = asarray([115.3,	230.8,	346.0,	461.5,	576.9,	691.2,	806.5],dtype='float64') * 1e9
nuij = Be * arange(1,nlevels) * 2
dvij = ones(len(nuij),dtype='float64') * 1e5  # delta-V = 1 km/s
# Einstein A values
mu = 0.1222e-18 # from http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html -> http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d028001.pdf -> http://adsabs.harvard.edu/abs/1994ApJS...95..535G
Aij = mu**2 * 64*pi**4 /(3*h*c**3) * nuij**4
#Aij = asarray([7.2e-8,	6.9e-7,	2.5e-6,	6.1e-6,	1.2e-5,	2.1e-5,	3.4e-5],dtype='float64')
#Aij = diag(A,1) #A is an NxN+1 array such that A[1,0]... just see the function A
# critical densities - not used in calcuations, but useful for comparison
ncr = asarray([1.1e3, 6.7e3, 2.1e4, 4.4e4, 7.8e4, 1.3e5, 2.0e5],dtype='float64')

def B_nu(frequency,T):
    """
    B_nu(T) - the Planck Function
    """
    T=float(T)
    myB_nu = 2 * h * frequency**3 / ( c**2 * ( exp( h*frequency / (k*T)) - 1 ) )
    return myB_nu

# not used Bij = diag(Bul,1) + diag(Blu,-1)
#Bij = Aij * c**2 / (2 * h * nuij**3)

#Cij is a function of density and temperature; density is solved for...
def Cul(n,T):
#    myCij = n * sigma * sqrt ( 2 * k * T / m_h2 )
    mygul = 4/sqrt(pi) * sqrt ( 2 * k * T / m_red ) * sigma * n         
    return mygul*ones(nuij.shape[0])

def Clu(n,T):
    myglu = Cul(n,T) * gi[1:]/gi[:-1] * exp( -1 * nuij * h / ( k * T ) )
    return myglu

def Gul(n,T):
    myGul = 4/sqrt(pi) * sqrt ( 2 * k * T / m_red ) * sigma * n
    return myGul * ones([nuij.shape[0]+1,nuij.shape[0]+1])
#* resize(fib,[fib.shape[0],nuij.shape[0]])

def Glu(n,T):
    myGlu = Gul(n,T) * append(gi[1:],gi[-1]+2)/gi * exp( -1 * nuij * h / ( k * T ) )
    return myGlu


# What is Rij?  Where's that written?
#assume beta, probability of escape, is 1
#def Rij(i,j,n,T):
#    if i < j:
#        myRij = Aij[i,j] * (1+Qij[i]) + Cij(n,T)
#    else:
#        myRij = gi[j]/gi[i] * Aij[i,j] * Qij[i] + Cij(n,T)
#    return myRij

#aRij = numpy.vectorize(Rij)

def A(n,tem,Rtem=2.73):
    """
    Matrix defining input/output from each level.  Includes collisions and radiation.
    """

    # allow for a higher-than-CMB radiation temperature (line trapping)
    Pij = B_nu(nuij,Rtem) # Planck function
    Q = c**2 / (2 * h * nuij**3) * Pij # radiation field at each point
    #this is not the einstein B coefficient: it's the einstein B coefficient times the rad field
    B = Aij * Q
    Bul = B 
    Blu = B * gi[1:] / gi[:-1]

    T=tem*ones(Aij.shape)
    size = Gul(n,tem).shape[0]
#    Amat_diag = diag( ones(Aij.shape[0]+1) + append(0,Cij(n,T)) + append(Cij(n,T),0) + append(0,Aij) +append(Blu,0) + append(0,Bul) )
#    Amat_diag = diag( append(0,Cij(n,T)) + append(Cij(n,T),0) + append(0,Aij) +append(Blu,0) + append(0,Bul) )
#modified: no collisional term for populating the zeroth state
    Amat_diag = diag( append(0,Cul(n,T)) + append(Clu(n,T),0) + append(0,Aij) +append(Blu,0) + append(0,Bul) )
#    Amat_diag =  diag( append(0,Aij) + append(Blu,0) + append(0,Bul) ) + identity(Gul(n,tem).shape[0],dtype='float64')*(Glu(n,tem)+Gul(n,tem))
#2/28    Amat_diag =  diag( append(0,Aij) + append(Blu,0) + append(0,Bul) ) + tri(size)*Gul(n,tem) + (1-tri(size))*Glu(n,tem)
#    Amat = diag(Bul,1) + diag(Blu,-1) + diag(Aij,1) + diag(Cij(n,T),1) + diag(Cij(n,T),-1) - Amat_diag
    Amat = diag(Bul,1) + diag(Blu,-1) + diag(Aij,1) + diag(Cul(n,T),1) + diag(Clu(n,T),-1) - Amat_diag
#    Amat = ( Gul(n,tem) + Glu(n,tem) ) * ( identity(Glu(n,tem).shape[0],dtype='float64')*-2+1) + diag(Bul,1) + diag(Blu,-1) + diag(Aij,1) - Amat_diag
#    Amat = ( Gul(n,tem) ) * ( identity(Glu(n,tem).shape[0],dtype='float64')*-1.+1.)  + diag(Bul,1) + diag(Blu,-1) + diag(Aij,1) - Amat_diag
#2/28    Amat = tri(size)*Glu(n,tem) + (1-tri(size))*Gul(n,tem)  + diag(Bul,1) + diag(Blu,-1) + diag(Aij,1) - Amat_diag
    newmat = zeros( asarray(Amat.shape) + 1 , dtype='float64')
#    newmat = zeros( [Amat.shape[0]+1,Amat.shape[1]], dtype='float64')
    newmat[:Amat.shape[0],:Amat.shape[1]] = Amat
    newmat[Amat.shape[0]] = 1
    return newmat

def ni(num,tem):
    """ num: density cm^-3 """
    b = zeros( Aij.shape[0] + 2 , dtype='float64')
    b[-1] = 10
    # this line is REALLY STUPID, but for some pointless reason linalg experiences 
    # precision errors and thinks that A is singular when it is obviously not.
    if tem > 30:
        myni = dot(linalg.pinv(A(num,tem)),b)
    else:
        scale = 1e6
        myni = linalg.solve( scale * (A(num,tem)) , b.T) / scale
#    myni = linalg.solve((A(num,tem)[:-1,:-1]), ones(Aij.shape[0]+1))
    return abs(myni)[:-1]

# T_ex = h nu / k  /  ( h nu / k T_kin + ln(1 + Aul Beta / Cul) )
def Tex(num,tem,trans,N=1e13):
    num = float64(num); tem = float64(tem)
    myni = abs(ni(num,tem))
#    myTex = h * nuij[trans] / k / log(myni[trans] * gi[trans+1] / (myni[trans+1] * gi[trans]) )
    myTex = h * nuij[trans] / k / ( h * nuij[trans] / (k * tem) + log( 1 + Aij[trans] * beta_trans(num,tem,N,trans) / Cul(num,tem)[trans] ) )
    return myTex

# Beta = ( 1 - e^-tau ) / tau   SPHERE
# Beta = ( 1 - e^-3 tau ) / 3 tau   PLANE / SLAB
# tau = Aul c^3 N_tot / ( 8 pi nu^3 delta-v )   * [ g_u/g_l * f_l + f_u ]
# f_u = n_u / n_tot
# delta-v = linewidth
# N_tot = column density

#niVect = numpy.vectorize(ni)
TexVect = numpy.vectorize(Tex)

def nx(tem):
    nunl = gi[1:]/gi[:-1] * exp(-h*nuij/ ( k * tem ) )
    return nunl

def tau(num,tem,N):
    N = float64(N)
    numerator = Aij * c**3 * N
    denominator = 8 * pi * nuij**3 * dvij
    myni = ni(num,tem)
    mynitot = myni.sum()
    factor = gi[1:]/gi[:-1] * myni[:-1] / mynitot + myni[1:] / mynitot
    return numerator / denominator * factor

def tau_trans(num,tem,N,trans):
    return tau(num,tem,N)[trans]

def beta(num,tem,N,shape='sphere'):
    mytau = tau(num,tem,N)
    if shape == 'sphere':
        return ( 1-exp(-mytau) ) / mytau
    elif shape == 'slab' or shape == 'plane':
        return ( 1-exp(-3*mytau) ) / ( 3 * mytau )

def beta_trans(num,tem,N,trans,shape='sphere'):
    return beta(num,tem,N,shape)[trans]

TauVect = numpy.vectorize(tau_trans)
BetaVect = numpy.vectorize(beta_trans)

def OM(T,freq=345e9):
    x0 =  k * T / ( h * freq ) 
    term2 = T * ( 1 - exp(-h*freq/(k*T)) ) 
    J = h * freq / k * (exp(h*freq/(k*T)) -1)**-1 
    return x0*term2*J


if __name__ == "__main__":

    x=logspace(0,6,1000)

    figure(8); clf();
    plot(logspace(12,17),[tau_trans(1e3,20,nn,0) for nn in logspace(12,17)])
    xlabel("Column Density (cm$^{-2}$)")
    ylabel("$\\tau$")
    savefig('/Users/adam/work/co/TauVsN.png',bbox_inches='tight')

    figure(9); clf()
    xlabel("Density (cm^-3)"); ylabel("Excitation Temperature (K)")
    title("20K gas")
    semilogx(x,TexVect(x,20,0),color='b',label="CO 1-0")
    semilogx(x,TexVect(x,20,1),color='g',label="CO 2-1")
    semilogx(x,TexVect(x,20,2),color='r',label="CO 3-2")
    semilogx(ncr[:3],TexVect(ncr,20,arange(7))[:3],'kx')
    gca().set_xlim(1,1e6)
    legend(loc='best')
    savefig('/Users/adam/work/co/CriticalDensity.png',bbox_inches='tight')

    figure(10); clf()
    xlabel("Density (cm^-3)"); ylabel("Excitation Temperature (K)")
    title("20K gas")
    semilogx(x,[Tex(i,20,0,N=1e12) for i in x],color='r',linestyle='--',label="N=1e12")
    semilogx(x,[Tex(i,20,0,N=1e13) for i in x],color='m',label="N=1e13")
    semilogx(x,[Tex(i,20,0,N=1e14) for i in x],color='k',label="N=1e14")
    semilogx(x,[Tex(i,20,0,N=1e15) for i in x],color='b',label="N=1e15")
    semilogx(x,[Tex(i,20,0,N=1e16) for i in x],color='g',label="N=1e16")
    semilogx(x,[Tex(i,20,0,N=1e17) for i in x],color='r',label="N=1e17")
    savefig('/Users/adam/work/co/CriticalDensityVsN.png',bbox_inches='tight')

    show()


