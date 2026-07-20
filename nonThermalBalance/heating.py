import numpy as np 
import astropy.units as u 

def heatingKasenBarnes(time_exp_days,averageAtomMass = 140,powerLaw = -1.3):
    '''
    Radioactive Power of r-process material.
    Power law of:
    Q = 1e10 t_d^-1.3  erg/s/g,
    where td is the time since merger measured in days from:
        https://iopscience.iop.org/article/10.3847/1538-4357/ab06c2/pdf  
    and references within.
    This function converts this to the work per ion, assuming an average atomic mass of  
    averageAtomMass nuclear units (140 by default).
    
    Using <A> = 140 with t = 1d, , this returns a work per ion of ~ 1.45 eV/s/ion,
    which is 10 * the number that Luke Shingles gave me.. - but need to multiply
    by a thermalization efficiency.
    '''

    const = 1e10 * u.erg / (u.g * u.s) #from 
    
    averageAtomicMass = (averageAtomMass * u.u).to('g')
    
    heating  = (averageAtomicMass * const).to('eV/s')
    
    return heating * (time_exp_days ** powerLaw)

def thermalizationKasenBarnes(time_exp_days,
                              ejectaMassSolar = 0.05,
                              eta = 1.0,
                              velocity_max_c=0.1,
                              temporalindex = 1.5,
                              fastelectron_parition = 0.2):
    
    ejectaMassSolar = ejectaMassSolar
    
    t_e = 15.0 * (eta * ejectaMassSolar / 0.01)**0.667 * (velocity_max_c / 0.2) **-2
    
    return fastelectron_parition * (1.0 + time_exp_days/t_e) ** (-1. * temporalindex)