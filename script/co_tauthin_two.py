import numpy as np
from astropy import units as u
from astropy import constants
from astroquery.splatalogue import Splatalogue

tcmb = 2.7315*u.K

# from https://spec.jpl.nasa.gov/ftp/pub/catalog/catdir.html
mu12co = 1.1011e-19 * u.esu*u.cm
mu13co = 1.1046e-19 * u.esu*u.cm
muc18o = 1.1079e-19 * u.esu*u.cm
B0_12CO = 57635.968 * u.MHz
B0_13CO = 55101.011 * u.MHz
B0_C18O = 54891.420 * u.MHz

# a way to get E_u
#Eu = Splatalogue.query_lines(100*u.GHz, 235*u.GHz, chemical_name=' CO ')['E_U (K)'][-1] * u.K

def J(T, nu=230.538*u.GHz):
    return (constants.h*nu/constants.k_B /
            (np.exp(constants.h * nu / (constants.k_B * u.Quantity(T, u.K))) -
             1)).to(u.K)

def prefactor(mu=mu12co, Ju=2, nu=230.538*u.GHz, B0=B0_12CO):
    return ((3 * constants.h) / (8 * np.pi**3 * mu**2 * Ju) * (constants.k_B / (constants.h * B0))).to((u.km/u.s)**-1 * u.cm**-2 * u.K**-1).value

def Nco_generic(Tex, flux, eu_k, nu, **kwargs):
    # comes from Mangum+ 2015 via https://github.com/keflavich/EtaCar_2013.1.00661.S/blob/master/analysis/co_xfactor.py
    Tex = u.Quantity(Tex, u.K)
    flux = u.Quantity(flux, u.K*u.km/u.s)
    return (prefactor(nu=nu, **kwargs) * (Tex.value+0.92)*np.exp(eu_k/Tex.value) *
            (np.exp(constants.h * nu / (constants.k_B * Tex)) - 1)**-1 * flux /
            (J(Tex) - J(tcmb)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2)

def Nco_13co_21(Tex, flux, eu_k=15.86632, nu=220.3986842*u.GHz, B0=B0_13CO, Ju=2, mu=mu13co, **kwargs):
    return Nco_generic(Tex=Tex, flux=flux, eu_k=eu_k, nu=nu, B0=B0, Ju=Ju, mu=mu, **kwargs) 

def Nco_12co_21(Tex, flux, eu_k=16.59608, nu=230.538*u.GHz, **kwargs):
    return Nco_generic(Tex=Tex, flux=flux, eu_k=eu_k, nu=nu, **kwargs) 
