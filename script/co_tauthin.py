"""
This script was started sometime around 2007 as part of a grad school homework
assignment, but I've expanded and reused it over time.  It still works,
and it makes some neat plots, but it is not an example of good coding practices.
This should all be rewritten with astropy units.
"""
from numpy import pi,exp,sqrt,log
import numpy

# CO 3-2 constants
hplanck = h = 6.626068e-27 # erg s 
k = kb = kB = 1.3806503e-16  # erg / K 
mu = 0.122e-18  # statcoulomb cm (.122 Debye, from NIST http://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Units=SI&Mask=1000#Diatomic ref 85)
               # also from muenther 1975
mu = 0.11011e-18 # from http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html -> http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d028001.pdf -> http://adsabs.harvard.edu/abs/1994ApJS...95..535G
# the above is mu_0, the rotationless dipole moment, which is not useful
mu = 0.1222e-18 # from http://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html -> http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d028001.pdf -> http://adsabs.harvard.edu/abs/1994ApJS...95..535G
coprops12CO10 = {'Aul':7.203e-8, 'Be': 57.63596828e9, 'nu': 115.271202e9 , 'Ju': 1}
coprops13CO10 = {'Aul':6.294e-8, 'Be': 55.1010138e9,  'nu': 110.201370e9 , 'Ju': 1}
coprops12CO21 = {'Aul':6.910e-7, 'Be': 57.63596828e9, 'nu': 230.5380e9,    'Ju': 2}
coprops13CO21 = {'Aul':6.038e-7, 'Be': 55.1010138e9,  'nu': 220.3986765e9, 'Ju': 2}
coprops12CO32 = {'Aul':2.497e-6, 'Be': 57.63596828e9, 'nu': 345.7959899e9, 'Ju': 3}
coprops13CO32 = {'Aul':2.181e-6, 'Be': 55.1010138e9,  'nu': 330.5879601e9, 'Ju': 3}
nu10 = 115.271202e9
nu32 = 345.7959899e9
msun = 1.99e33 # grams
level3pop = .176  # from IP1 spreadsheet;  not exp(-h * 345e9 / (k * 20))
h2toCO    = 1e4 # X_CO = CO abundance
kmstocms  = 1e5 # can't forget to get EVERYTHING in CGS
eta_mb    = .68 # Main Beam Efficiency for HARP-B on JCMT using Saturn
eta_mb    = .60 # a "typical" main beam efficiency for Mars http://www.jach.hawaii.edu/JCMT/spectral_line/Standards/current_cals.html
vcen      = -16.0 # approximate velocity of cloud
c = 2.99792458e10  # cm/s
Be=57.64e9     # CO 3-2 (NIST)
Aul = A32 = 2.5e-6   # CO 3-2 (Turner 1977)
Aul = A32 = 2.497e-06   # CO 3-2 (Leiden molecular database: http://www.strw.leidenuniv.nl/~moldata/datafiles/co.dat)
Jl32 = 2
A32_mu = 64 * pi**4 / (3*h*c**3) * nu32**3 * (Jl32+1.0)/(2*Jl32+3) * mu**2 # Rohlfs and Wilson 15.24
A10 = 7.203e-8 # (Leiden molecular database: http://www.strw.leidenuniv.nl/~moldata/datafiles/co.dat)
Jl10 = 0
A10_mu = 64 * pi**4 / (3*h*c**3) * nu10**3 * (Jl10+1.0)/(2*Jl10+3) * mu**2 # Rohlfs and Wilson 15.24
#print "A10_mu / A10" , A10_mu/A10
#print "A32_mu / A32" , A32_mu/A32
mh=1.67262158e-24
mumass = 1.4
massmu = 1.4
mh2 = 2*mh
Tex = 20
tcmb = 2.73
dist_kpc = 2.0
distpc = dist_pc = dist_kpc*1.0e3  # distance to W5 
distance = dist_pc
au = 1.496e13
pc = 3.086e18
muh2=2.8
yr = year = 86400*365.

Tex = 20
mmtocol = 2.19e22 * (exp(13.01/Tex)-1)
mmtoav = mmtocol / 1.9e20
mmtomass = mmtocol * muh2 * mh / msun # 14.30 * (exp(13.01/Tex) - 1) * dist_kpc**2
mmtomass = 14.30 * (exp(13.01/Tex) - 1) * dist_kpc**2


fwhm_const = numpy.sqrt(8*numpy.log(2))
fwhm_bgps =33.0 # arcsec
fwhm_fcrao=46.0 # arcsec
omegabeam_bgps  = pi*((fwhm_bgps /fwhm_const*sqrt(2))/206265.0)**2
omegabeam_fcrao = pi*((fwhm_fcrao/fwhm_const*sqrt(2))/206265.0)**2
ppbeam_bgps =  pi * ( sqrt(2)/fwhm_const * fwhm_bgps  )**2 / 7.2**2  # pi r_eff^2 / x^2
ppbeam_fcrao = pi * ( sqrt(2)/fwhm_const * fwhm_fcrao )**2 / 7.2**2 


def nj(tem,Ju,numlevs=50,Be=Be):
    J = numpy.arange(numlevs)
    ntotDn0 = sum( (2*J+1)*exp(-J*(J+1.0)*Be*h/(kb*tem)) )
    #ntotDn0 = sum( (2*J+1)*exp(-J*(J+1.0)*Be*h/(kb*tem)) - (2*J+1)*exp(-J*(J+1.0)*Be*h/(kb*tcmb)) )
    if Ju == 0:
        return 1.0/ntotDn0
    expBe = exp(-h*Be*(Ju*(Ju+1))/ ( kb * tem ) )
    nuDn0 = (2*Ju+1.0) * expBe
    #ntotDnu = ntotDn0 * (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( kb * tem ) ) 
    #ntotDnu = ntotDn0 * ((2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( kb * tem ) ) - \
    #                    (2.0*Jl+1.0)/(2.0*Ju+1.0) * exp(h*Be*(Ju*(Ju+1))/ ( kb * tcmb ) ) )
    return nuDn0/ntotDn0

def njapprox(tem,Ju,Be=Be):
    ntotDn0 = kb * tem / (Be*h)
    # ntot/n0 = ntotDn0... n0 = ntot/ntotDn0
    # nj/n0 = gj/g0 * exp(-Ej0 / kB T)
    # ntot / nj = ntotDn0 / (gj/j0) * exp(Ej0/kbT)
    if Ju == 0:
        return 1.0/ntotDn0
    ntotDnu = ntotDn0 / ( (2.0*Ju+1.0) * exp(-h*Be*(Ju*(Ju+1))/ ( kb * tem ) ) ) 
    return 1.0/ntotDnu

def cotocolapprox(tex,coprops,conversionratio=h2toCO,old=False):
    """
    Conversion ratio - ratio of particular CO molecule to H2
    """
    nu = coprops['nu']
    Ju = coprops['Ju']
    Aul = coprops['Aul']
    Be = coprops['Be']
    if not isinstance(tex,numpy.ndarray):
        tex = numpy.array(tex)
    if old:
        return conversionratio * (8*pi*nu**3*kb)/((2*Ju+1)*c**3*h*Be*Aul) * exp(Ju*(Ju+1)*h*Be/(kb*tex)) * (exp(h*nu/(kb*tex))-1)**-1 * kmstocms # From eqn 1 of paper
    else:
        ntotDnu = 1.0 / numpy.array([ njapprox(T,Ju,Be=Be) for T in tex ]) 
        #numupper = (8*pi*nu**3)/(c**3*Aul*tex) * (exp(h*nu/(kb*tex))-1)**-1 * kmstocms # From eqn 1 of paper
        expcmb = exp(h*nu/(kb*tcmb))
        numupper = ( 8 * pi * nu**2 ) / ( c**3 * Aul ) * kb / h * (expcmb - 1) / ( expcmb - (exp(hplanck*nu/(kb*tex))) ) # - (exp(hplanck*nu/(kb*tcmb))-1)**-1)
        return conversionratio * ntotDnu * numupper * kmstocms

def cotocol(tex,coprops,conversionratio=h2toCO):
    """
    Conversion ratio - ratio of particular CO molecule to H2
    """
    nu = coprops['nu']
    Ju = coprops['Ju']
    Be = coprops['Be']
    Aul = coprops['Aul']
    if not isinstance(tex,numpy.ndarray):
        tex = numpy.array(tex)
    ntotDnu = 1.0 / numpy.array([ nj(T,Ju,Be=Be,numlevs=50) for T in tex ]) 
    expcmb = exp(h*nu/(kb*tcmb))
    Nupper = ( 8 * pi * nu**2 ) / ( c**3 * Aul ) * kb / h * (expcmb - 1) / ( expcmb - (exp(hplanck*nu/(kb*tex))) ) # - (exp(hplanck*nu/(kb*tcmb))-1)**-1)
    return conversionratio * ntotDnu * Nupper * kmstocms # From eqn 1 of paper

def cotocol_nocmb(tex,coprops,conversionratio=h2toCO):
    """
    Conversion ratio - ratio of particular CO molecule to H2
    """
    nu = coprops['nu']
    Ju = coprops['Ju']
    Be = coprops['Be']
    Aul = coprops['Aul']
    if not isinstance(tex,numpy.ndarray):
        tex = numpy.array(tex)
    ntotDnu = 1.0 / numpy.array([ nj(T,Ju,Be=Be,numlevs=50) for T in tex ]) 
    Nupper = ( 8 * pi * nu**2 ) / ( c**3 * Aul ) * kb / h 
    return conversionratio * ntotDnu * Nupper * kmstocms # From eqn 1 of paper


def taucotocol(tbeam,tex,coprops,conversionratio=h2toCO):
    """
    tbeam - measured antenna temperature
    Conversion ratio - ratio of particular CO molecule to H2
    
    (uses eqn 19 of column_derivation (as of Oct 5, 2010)) 
    """
    nu = coprops['nu']
    Ju = coprops['Ju']
    Be = coprops['Be']
    Aul = coprops['Aul']
    ntotDnu = 1.0 / numpy.array([ nj(T,Ju,Be=Be,numlevs=50) for T in tex ]) 
    expcmb = exp(h*nu/(kb*tcmb))
    exptex = exp(h*nu/(kb*tex))
    tau = -log(1-kb*tbeam/(h*nu) * ( (exptex - 1)**-1 - (expcmb-1)**-1 )**-1 )
    Nupper = ( 8 * pi * nu**3 ) / ( c**3 * Aul ) * (exptex-1)**-1 * tau
    return conversionratio * ntotDnu * Nupper * kmstocms # From eqn 1 of paper

cotocol13co   = 1.694e20 / (1-exp(-5.3/Tex))  # integrated 13co to column density conversion 
cotocol13co10 = cotocol(numpy.array([20]),coprops13CO10,conversionratio=6e5)
cotocol32     = cotocol(numpy.array([20]),coprops12CO32) #h2toCO * (8*pi*nu32**3*kb)/(7*c**3*h*Be*Aul) * exp(4*3*h*Be/(kb*Tex)) * (exp(h*nu32/(kb*Tex))-1)**-1 * kmstocms # From eqn 1 of paper
cotocol10     = cotocol(numpy.array([20]),coprops12CO10) #h2toCO * (8*pi*nu10**3*kb)/(3*c**3*h*Be*A10) * exp(2*1*h*Be/(kb*Tex)) * (exp(h*nu10/(kb*Tex))-1)**-1 * kmstocms # From eqn 1 of paper
# also, constants @ http://casa.colorado.edu/~ginsbura/COproperties.htm
beamarea32_cm2 = 2*pi*(18.0/2.35 * dist_pc * au)**2
cotomass32 = cotocol32 * massmu * mh2 / msun * beamarea32_cm2
beamarea10_cm2 = 2*pi*(45.0/2.35 * dist_pc * au)**2
cotomass13co10 = cotocol13co10 * massmu * mh2 / msun * beamarea10_cm2
# for T_ex = 20K, cotocol = 7.28e20
# cotomass = 0.275 msun/pixel for 22.5" square pixels

if __name__=="__main__":
    from pylab import *
    print("Plotting up CO vs. Tex for debug purposes")
    texarr = linspace(3,60,200)
    twelveCO10   = cotocol(texarr,coprops12CO10)
    twelveCO21   = cotocol(texarr,coprops12CO21)
    twelveCO32   = cotocol(texarr,coprops12CO32)
    thirteenCO10 = cotocol(texarr,coprops13CO10,conversionratio=6e5)
    thirteenCO21 = cotocol(texarr,coprops13CO21,conversionratio=6e5)
    thirteenCO32 = cotocol(texarr,coprops13CO32,conversionratio=6e5)
    twelveCO10approx   = cotocolapprox(texarr,coprops12CO10)
    twelveCO21approx   = cotocolapprox(texarr,coprops12CO21)
    twelveCO32approx   = cotocolapprox(texarr,coprops12CO32)
    thirteenCO10approx = cotocolapprox(texarr,coprops13CO10,conversionratio=6e5)
    twelveCO10approxold   = cotocolapprox(texarr,coprops12CO10,old=True)
    twelveCO32approxold   = cotocolapprox(texarr,coprops12CO32,old=True)
    thirteenCO10approxold = cotocolapprox(texarr,coprops13CO10,old=True,conversionratio=6e5)
    twelveCO10_nocmb   = cotocol_nocmb(texarr,coprops12CO10)
    twelveCO21_nocmb   = cotocol_nocmb(texarr,coprops12CO21)
    twelveCO32_nocmb   = cotocol_nocmb(texarr,coprops12CO32)

    # Arce et al 2010 and Arce and Goodman 2001
    T0_13 = h*coprops13CO10['nu']/kb
    T0_12 = h*coprops12CO10['nu']/kb
    hb3k_13 = h*coprops13CO10['Be']/(3*kb)
    hb3k_12 = h*coprops12CO10['Be']/(3*kb)
    #arce13co = 2.42e14*(texarr+hb3k_13) / (1-exp(-h*coprops13CO10['nu']/(kb*texarr))) /(T0_13/(exp(T0_13/texarr)-1)-T0_13/(exp(T0_13/tcmb)-1)) * 7e5
    arce13co = ((3*h/(8*pi**3*mu**2)*kb*kmstocms/(h*coprops13CO10['Be'])))*(texarr+hb3k_13) / (1-exp(-h*coprops13CO10['nu']/(kb*texarr))) /(T0_13/(exp(T0_13/texarr)-1)-T0_13/(exp(T0_13/tcmb)-1)) * 7e5
    arce12co = ((3*h/(8*pi**3*mu**2)*kb*kmstocms/(h*coprops12CO10['Be'])))*(texarr+hb3k_12) / (1-exp(-h*coprops12CO10['nu']/(kb*texarr))) /(T0_12/(exp(T0_12/texarr)-1)-T0_12/(exp(T0_12/tcmb)-1)) * h2toCO
    # Cabrit and Bertout 1990
    # Where does their "Omega" factor come from?  It looks a lot like the radiative transfer
    # equation but it is not...
    nu10 = coprops12CO10['nu']
    n0_cabrit = 1.0/(kb * texarr / (h * Be))
    omega_cabrit = (n0_cabrit * (1-exp(-h*nu10/(kb*texarr))) * \
            ( (h*nu10/kb) * ( (exp(h*nu10/(kb*texarr))-1)**-1 - (exp(h*nu10/(kb*tcmb))-1)**-1) ) )
    cabrit12co = h2toCO * 3 * h / (8 * pi**3 * (mu)**2) / omega_cabrit * kmstocms

    # lada & fich 1996
    n0_lada_11 = 1.0/array([ nj(T,0,numlevs=11) for T in texarr ])
    omega_lada = (n0_lada_11 * (1-exp(-h*nu10/(kb*texarr))) * \
            ( h*nu10/kb * ( (exp(h*nu10/(kb*texarr))-1)**-1 - (exp(h*nu10/(kb*tcmb))-1)**-1) ) )

    n0_lada_100 = 1.0/array([ nj(T,0,numlevs=100) for T in texarr ])
    omega_lada_100 = (n0_lada_100 * (1-exp(-h*nu10/(kb*texarr))) * \
            ( h*nu10/kb * ( (exp(h*nu10/(kb*texarr))-1)**-1 - (exp(h*nu10/(kb*tcmb))-1)**-1) ) )
    lada12co100 =  h2toCO * 3 * h / (8 * pi**3 * (mu)**2) / omega_lada_100 * kmstocms 

    n1_8      = array([ nj(T,1,numlevs=8 ) for T in texarr ]) 
    n1_10     = array([ nj(T,1,numlevs=10) for T in texarr ]) 
    n1_50     = array([ nj(T,1,numlevs=50) for T in texarr ])
    n1_100    = array([ nj(T,1,numlevs=100) for T in texarr ])
    n1_approx = array([ njapprox(T,1) for T in texarr ])
    n0_100    = array([ nj(T,0,numlevs=100) for T in texarr ])
    n0_approx = array([ njapprox(T,0) for T in texarr ])

    figure(4); clf()
    plot(texarr,n0_100,label='100')
    plot(texarr,n0_approx,label='approx')
    legend(loc='best')
    xlabel('Excitation Temperature')
    ylabel('$n_{tot}/n_0$')
    savefig('tex_vs_density.png', bbox_inches='tight')

    figure(3); clf();
    plot(texarr,n1_8     ,'--',label='8',linewidth=3)
    plot(texarr,n1_10    ,'--',label='10',linewidth=2.5)
    plot(texarr,n1_50    ,'--',label='50',linewidth=2)
    plot(texarr,n1_100   ,'--',label='100',linewidth=1.5)
    plot(texarr,n1_approx,'--',label='approx')
    plot(texarr,omega_cabrit,label='cabrit')
    plot(texarr,omega_lada,label='lada')
    plot(texarr,omega_lada_100,label='lada_100')
    legend(loc='best')
    xlabel('Excitation Temperature')
    ylabel('$\\Omega$')
    savefig('tex_vs_Omega.png', bbox_inches='tight')

    figure(1); clf();
    plot1 = semilogy(texarr,twelveCO10,label='Exact partition function',linewidth=2.0)
    ax = plot1[0].axes
    semilogy(texarr,twelveCO10approx,label='Approximate paritition function',linewidth=2.0)
    legend(loc='best')
    xlabel('Excitation Temperature')
    ylabel('N(H$_2$) / K km s$^{-1}$')
    ax.set_ylim(8e18,3e19)
    savefig('/Users/adam/work/co/column_derivation/columnconversion_approximation.png',bbox_inches='tight')
    savefig('columnconversion_approximation.png',bbox_inches='tight')
    #semilogy(texarr,twelveCO32,'--',label='12CO32')
    #semilogy(texarr,twelveCO32approx,':',label='12CO32 approx')
    #semilogy(texarr,twelveCO32approxold,'-',label='12CO32 approxold')
    plot(texarr,cabrit12co,'k:',label='cabrit12co',linewidth=2.0)
    plot(texarr,lada12co100,'c:',label='lada12co100',linewidth=2.0)
    plot(texarr,arce12co,'--',label='arce12co',linewidth=2.0)
    #plot(texarr,arce13co/70.,'-.',label='arce13co',linewidth=2.0)
    #semilogy(texarr,twelveCO10approxold,'--',label='12CO10 approxold')
    legend(loc='best')
    ax.set_ylim(8e18,3e19)
    savefig('columnconversion_approximation_withothers.png', bbox_inches='tight')

    figure(2); clf();
    tbeam = linspace(0.01,50,100)
    col_12co_tau_20K = taucotocol(tbeam,numpy.array([20]),coprops12CO10)
    col_12co_tau_50K = taucotocol(tbeam,numpy.array([50]),coprops12CO10)
    col_12co_approx_20K = tbeam*cotocol(numpy.array([20]),coprops12CO10)
    col_12co_approx_50K = tbeam*cotocol(numpy.array([50]),coprops12CO10)
    plot(tbeam,col_12co_tau_20K/col_12co_approx_20K,label='T$_{ex}$=20K')
    plot(tbeam,col_12co_tau_50K/col_12co_approx_50K,label='T$_{ex}$=50K')
    #semilogy(tbeam,col_12co_approx_20K,label='20K approx')
    #semilogy(tbeam,col_12co_approx_50K,label='50K approx')
    xlabel('Beam Temperature')
    ylabel('Ratio N(H$_2$): exact / approximate')
    legend(loc='best')
    grid()
    gca().set_ylim(1,2)
    savefig('/Users/adam/work/co/column_derivation/ratio_exactapprox.png',bbox_inches='tight')
    savefig('ratio_exactapprox.png',bbox_inches='tight')

    figure(5); clf();
    plot2 = semilogy(texarr,twelveCO32,'-',label='12CO32',linewidth=3.0)
    ax2 = plot2[0].axes
    semilogy(texarr,twelveCO32approx,'--',label='12CO32 approx',linewidth=2.0)
    ax2.set_ylim(3e18,5e19)
    legend(loc='best')
    xlabel('Excitation Temperature')
    ylabel('N(H$_2$) / K km s$^{-1}$')
    savefig('/Users/adam/work/co/column_derivation/co32conversion.png',bbox_inches='tight')
    scatter(50,2.5e19,label='Hatchell 2007')
    ax2.set_ylim(3e18,5e19)
    savefig('/Users/adam/work/co/column_derivation/co32conversion_hatchellcompare.png',bbox_inches='tight')
    savefig('co32conversion_hatchellcompare.png',bbox_inches='tight')

    rc("font",size=24)
    figure(6); clf();
    plot3 = semilogy(texarr,twelveCO10,label='$^{12}$CO 1-0',linewidth=2.0,color='b')
    ax3 = plot3[0].axes                                                 
    semilogy(texarr,twelveCO21,label='$^{12}$CO 2-1',linewidth=2.0,color='g')
    semilogy(texarr,twelveCO32,label='$^{12}$CO 3-2',linewidth=2.0,color='r')
    leg=legend(loc='best')
    [T.set_fontsize(14) for T in leg.get_texts()]
    xlabel('Excitation Temperature (K)')
    ylabel('N(H$_2$) / K km s$^{-1}$')
    ax3.set_ylim(3e18,5e19)
    ax3.set_yticks([4e18,1e19,3e19])
    ax3.set_yticklabels(["$%i\\times10^{%i}$" % (x,y) for x,y in [(4,18),(1,19),(3,19)]])
    savefig('/Users/adam/work/co/column_derivation/12CO_columnconversion_vs_tex.png',bbox_inches='tight')
    savefig('12CO_columnconversion_vs_tex.png',bbox_inches='tight')

    #semilogy(texarr,thirteenCO10/60.,label='$^{13}$CO 1-0',linewidth=2.0,linestyle='dashed',color='b')
    #semilogy(texarr,thirteenCO21/60.,label='$^{13}$CO 2-1',linewidth=2.0,linestyle='dashed',color='g')
    #semilogy(texarr,thirteenCO32/60.,label='$^{13}$CO 3-2',linewidth=2.0,linestyle='dashed',color='r')
    semilogy(texarr,twelveCO10approx,label='CO 1-0 approx',linewidth=2.0,linestyle='dashed',color='b')
    semilogy(texarr,twelveCO21approx,label='CO 2-1 approx',linewidth=2.0,linestyle='dashed',color='g')
    semilogy(texarr,twelveCO32approx,label='CO 3-2 approx',linewidth=2.0,linestyle='dashed',color='r')
    #legend(loc='best')
    xlabel('Excitation Temperature (K)')
    ylabel('N(H$_2$) / K km s$^{-1}$')
    ax3.set_ylim(3e18,5e19)
    ax3.set_yticks([4e18,8e18,1e19,2e19,4e19])
    ax3.set_yticklabels(["$"+("%i\\times" % x)*z+"10^{%i}$" % y for x,y,z in [(4,18,1),(8,18,1),(1,19,0),(2,19,1),(4,19,1)]])
    draw()
    savefig('/Users/adam/work/co/column_derivation/columnconversion_vs_tex.png',bbox_inches='tight')
    savefig('columnconversion_vs_tex.png',bbox_inches='tight')

    semilogy(texarr,twelveCO10_nocmb,'b:',linewidth=2.0)
    semilogy(texarr,twelveCO21_nocmb,'g:',linewidth=2.0)
    semilogy(texarr,twelveCO32_nocmb,'r:',linewidth=2.0)
    ax3.set_ylim(3e18,5e19)
    leg=legend(loc='best')
    savefig('/Users/adam/work/co/column_derivation/columnconversion_vs_tex_allapprox_legend.png',bbox_inches='tight')
    savefig('columnconversion_vs_tex_allapprox_legend.png',bbox_inches='tight')
    leg.set_visible(False)
    savefig('/Users/adam/work/co/column_derivation/columnconversion_vs_tex_allapprox.png',bbox_inches='tight')
    savefig('columnconversion_vs_tex_allapprox.png',bbox_inches='tight')

    temperatures = arange(12)
    njcmb = array([nj(array([2.73]),Ju,200) for Ju in temperatures])
    nj5   = array([nj(array([5]   ),Ju,200) for Ju in temperatures])
    nj10  = array([nj(array([10]  ),Ju,200) for Ju in temperatures])
    nj15  = array([nj(array([15]  ),Ju,200) for Ju in temperatures])
    nj20  = array([nj(array([20]  ),Ju,200) for Ju in temperatures])
    nj30  = array([nj(array([30]  ),Ju,200) for Ju in temperatures])
    nj50  = array([nj(array([50]  ),Ju,200) for Ju in temperatures])
    nj100 = array([nj(array([100] ),Ju,200) for Ju in temperatures])

    figure(7); clf();
    plot(temperatures,njcmb, 's--',label='2.73K')
    plot(temperatures,nj5 , 's-',label='5K')
    plot(temperatures,nj10, 's-',label='10K')
    #plot(temperatures,nj15, 's-',label='15K')
    plot(temperatures,nj20, 's-',label='20K')
    plot(temperatures,nj30, 's-',label='30K')
    plot(temperatures,nj50, 's-',label='50K')
    #plot(temperatures,nj100,'s-',label='100K')
    legend(loc='best')
    xlabel('Upper rotational level of CO')
    ylabel('Fractional level population')
    savefig('/Users/adam/work/co/column_derivation/fractional_level_population.png',bbox_inches='tight')
    savefig('fractional_level_population.png',bbox_inches='tight')

    n0ofT = array([ nj(T,0,200) for T in texarr ])
    n1ofT = array([ nj(T,1,200) for T in texarr ])
    n2ofT = array([ nj(T,2,200) for T in texarr ])
    n3ofT = array([ nj(T,3,200) for T in texarr ])

    figure(8); clf();
    plot(texarr,n0ofT,label='J=0')
    plot(texarr,n1ofT,label='J=1')
    plot(texarr,n2ofT,label='J=2')
    plot(texarr,n3ofT,label='J=3')
    legend(loc='best')
    xlabel('Excitation Temperature')
    ylabel('Fractional level population')
    xticks(arange(0,60,5))
    grid()
    savefig('/Users/adam/work/co/column_derivation/fractional_level_population_vs_temp.png',bbox_inches='tight')
    savefig('fractional_level_population_vs_temp.png',bbox_inches='tight')

    figure(9); clf();
    sp1=subplot(211)
    semilogy(texarr,twelveCO10,'r',label='CO 1-0')
    semilogy(texarr,twelveCO21,'g',label='CO 2-1')
    semilogy(texarr,twelveCO32,'b',label='CO 3-2')
    semilogy(texarr,twelveCO10_nocmb,'r--')
    semilogy(texarr,twelveCO21_nocmb,'g--')
    semilogy(texarr,twelveCO32_nocmb,'b--')
    legend(loc='best')
    xlabel('T$_{ex}$ (K)')
    ylabel('N(H$_2$) [cm$^{-2}$] / (K km s$^{-1}$)')
    sp1.set_ylim(4e18,3e19)
    sp2=subplot(212)
    semilogy(texarr,log10(twelveCO10/twelveCO10_nocmb),'r',label='CO 1-0')
    semilogy(texarr,log10(twelveCO21/twelveCO21_nocmb),'g',label='CO 2-1')
    semilogy(texarr,log10(twelveCO32/twelveCO32_nocmb),'b',label='CO 3-2')
    xlabel('T$_{ex}$ (K)')
    ylabel('Correction Factor')
    sp2.set_xticks(arange(0,60,5))
    ticks = array([1.0002,1.0005,1.001,1.002,1.005,1.01,1.02,1.05])
    sp2.set_yticks(ticks-1) #logspace(-4,-0.7,11))
    sp2.set_yticklabels(ticks) #["%5.4f" % t for t in (1+logspace(-4,-0.7,11))])
    grid()
    sp2.set_ylim(0,log10(1.2))
    savefig('/Users/adam/work/co/column_derivation/CMB_correction_factor.png',bbox_inches='tight')
    savefig('CMB_correction_factor.png',bbox_inches='tight')

    draw()

    #show()

    
