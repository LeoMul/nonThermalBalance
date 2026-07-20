'''
Calculates the the probability P_{ij} that the photon from a  recombination event 

X^{j+1} + e^- ---> X^{j} + hv 

photoionizes X^{i}. This is importance as the prevalence of higher charge states
in the presence of non-thermal electrons can in principle further undermine the 
near-neutral populations.

'''


#imports

import numpy as np 
import astropy.units as u 
from constants import PHI_R_DEFAULT

def pijAxelrod(fractions, pi_thresholds, phi_r = PHI_R_DEFAULT , coldens = np.inf* u.cm**-2):
    
    #Dimension of self-photoionizations possible
    nst = len(fractions) - 1 
    
    pij = np.zeros([nst,nst])
    
    for j in range(0, nst):
        sumoverk = 1e-30 * u.cm * u.cm 
        
        for k in range(0,j+1): #python looping
            sumoverk += fractions[k] * pi_thresholds[k,j]
        
        tau_j = coldens * sumoverk
        tau_j = (1.0 - np.exp(-1. * tau_j))
        
        for i in range(0, nst):
            pij[i,j] = tau_j * fractions[i]*pi_thresholds[i,j] / sumoverk
            
    return pij * phi_r

def pijHotokezaka(fractions,temp, PICrossList, ionizationPotentials):
    
    from scipy.interpolate import interp1d
    
    nst = len(fractions) - 1 

    energyGrid = np.linspace(0,150,50000) * u.eV
    
    totalCross = np.zeros_like(energyGrid.value) * u.cm * u.cm
    
    
    partialCross = np.zeros([nst, len(energyGrid)])* u.cm * u.cm
    
    counter = 0 
    
    for cross in PICrossList:
        thisGrid = cross[:,0]
        thisCros = cross[:,1]
        select = energyGrid.value > thisGrid[0]
        
        interp = interp1d(thisGrid,thisCros)
        
        thisint = interp(energyGrid[select])*u.cm*u.cm
        
        totalCross[select] += thisint
        
        partialCross[counter][select] = thisint
        counter+= 1 
        

    
    pij = np.zeros([nst,nst])


    totalCrossFraction = np.zeros_like(energyGrid.value) * u.cm * u.cm 
    partialCrossDivided = (partialCross.value).copy()

    for ii in range(0,nst):
        totalCrossFraction = totalCrossFraction + fractions[ii] * partialCross[ii,:]

    for ii in range(0,nst):
        partialCrossDivided[ii,:] *= fractions[ii]
        partialCrossDivided[ii,:] /= (totalCrossFraction.value+1e-30)

    norm = 1. 


    for jj in range(0,nst):
        
        #right now normalized s.t recombination energy = I_j + kT 
        
        #norm = 1. + ionizationPotentials[jj] / (temp * c.k_B)

        select =  (energyGrid <= ionizationPotentials[jj] + temp * c.k_B) #& (energyGrid >= ionizationPotentials[jj])

        egridJ = energyGrid[select]

        dndhv_j = norm/ ( egridJ + 1e-30*u.eV)

        for ii in range(0,nst):
            pij[ii,jj] = np.trapezoid(dndhv_j * partialCrossDivided[ii][select], egridJ)
            
    return pij 




def getPhotoionizationCrossSections():
    
    
    return None 