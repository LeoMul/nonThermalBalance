
import numpy as np 
from photonRecycling import * 
from constants import *
MAX_PHOTO_ITER = 20 #max iter for self-photoionization.

def ionizationBalance(ionization_rates, recombination_rates,atomicNumber, photonRecycling=True,pi_thresholds=None, phi_r = PHI_R_DEFAULT, verbose = True):
    ''' 
    Coronal ionization balance. 
    I.e - only consider:
        - ionization    out of ground into all adjacent states
        - recombination out of ground into all adjacent states
    
    ionization_rates   [i] = rate of ionization    from i   to i+1 
    recombination_rates[i] = rate of recombination from i+1 to i    
    
    '''    

    populations = ionizationBalanceForwardIter(ionization_rates, recombination_rates,atomicNumber)
    
    #hack, for now
    
    PICrossList = []

    ll = ['/Users/leomulholland/CePaper/PIData/CeI/xout1','/Users/leomulholland/CePaper/PIData/CeII/xout1','/Users/leomulholland/CePaper/PIData/CeIII/xout1','/Users/leomulholland/CePaper/PIData/CeIV/xout1','/Users/leomulholland/CePaper/PIData/CeV/xout1']
    ionizationPotentials = np.zeros(5) * u.eV
    counter = 0
    for l in ll:
        PICross = np.loadtxt(l,usecols=(1,2))
        PICross[:,0] = (PICross[:,0] * u.rydberg).to("eV")
        ionizationPotentials[counter] = PICross[0,0] * u.eV
        PICross[:,1] = PICross[:,1] * 1e-18 * u.cm * u.cm
        PICrossList.append(PICross)
        counter += 1 
    pi_thresholds = np.zeros([5,5]) * u.cm * u.cm
    from scipy.interpolate import interp1d
    for ii in range(0,5):
        PICross = PICrossList[ii]
        thisCrossInterp = interp1d(PICross[:,0], PICross[:,1]) 
        for jj in range(ii,5):


            pi_thresholds[ii,jj] = thisCrossInterp(ionizationPotentials[jj])* u.cm * u.cm
    
    #print(pi_thresholds)
    
    
    
    if (photonRecycling and (pi_thresholds is None)):
        print("Warning, photonRecycling is on but I have no PI data. Setting all PI cross sections to 1 Mb = 1e-18 cm^2.")
        nst = len(populations)
        pi_thresholds = np.zeros([nst,nst]) * u.cm**2
        for ii in range(0,nst):
            for jj in range(ii,nst):
                pi_thresholds[ii,jj] = PI_CROSS_THRESH_CM2_DEFAULT * u.cm ** 2
    
        
    populationsNoPhoto  = populations.copy()
    
    if photonRecycling:
        for _ in range(0, MAX_PHOTO_ITER):
            pij = pijAxelrod(populations, pi_thresholds, phi_r)
            
            populations = ionizationBalanceForwardIter(ionization_rates, recombination_rates,atomicNumber,pij, populations)


    if verbose:
        nstg = len(populations)
        print(' Balance:')
        print('    NoPR    Pr')
        
        for ii in range(0,nstg):
            print('    {:7.2e}    {:7.2e} '.format(populationsNoPhoto[ii], populations[ii]))


    
    return populations


def ionizationBalanceForwardIter(ionization_rates, recombination_rates,atomicNumber, pij = None, f_prev = None):
    ''' 
    Coronal ionization balance. 
    I.e - only consider:
        - ionization    out of ground into all adjacent states
        - recombination out of ground into all adjacent states
    
    ionization_rates   [i] = rate of ionization    from i   to i+1 
    recombination_rates[i] = rate of recombination from i+1 to i    
    
    '''    
    
    N = len(ionization_rates)
    num_states = N + 1
    
    
    f_new = np.zeros(num_states)
    
    # Set the ground state (neutral fraction) as our relative baseline
    f_new[0] = 1.0
    

    for i in range(0,N):
        # Calculate the photoionization loss field from higher states
        photo_field = 0.0
        pii         = 0.0 
        if ( (pij is not None) and (f_prev is not None) ):
            #
            pii = pij[i, i]
            #
            for j in range(i + 1, N-1):
                photo_field += pij[i, j] * recombination_rates[j] * f_prev[j + 1]
            #    
            photo_field/= f_prev[i]
    
        numerator = (ionization_rates[i] + photo_field )  * f_new[i] 

        denominator = (1.0 - pii) * recombination_rates[i]
        
        #print(numerator, denominator)
        
        f_new[i + 1] = numerator / denominator
    
    #Top up with geomtric seris on the high charge states not included 
    #really this is just to check if the user is missing high charge states.
    pp = f_new[-1]
    norm = f_new.sum()
    normOld = norm.copy()
    
    current = pp 
    totalarray = [*f_new]
    for _ in range(len(ionization_rates),atomicNumber):
        current = current * 0.25 
        norm += current 
        totalarray.append(current)
    
    #print(normOld,norm)
    
    if (norm / normOld -1 > IONBALANCE_NORM_TOLERANCE):
        import sys
        print('Ionization balance norm is not converged, consider adding more ionization stages.') 
        sys.exit()
    
    
    totalarray = np.array(totalarray)
    totalarray /= norm
    f_new = f_new / norm
    
    return f_new



#Retired Method:

def ionizationBalanceOld(ionization_rates, recombination_rates,atomicNumber):
    ''' 
    Coronal ionization balance. 
    I.e - only consider:
        - ionization    out of ground into all adjacent states
        - recombination out of ground into all adjacent states
    
    ionization_rates   [i] = rate of ionization    from i   to i+1 
    recombination_rates[i] = rate of recombination from i+1 to i    
    
    pads out the balance, by assuming the ionization decreases 
    geometrically from the last explicitly calculated and that the 
    recombination increases geometrically. 
    
    assuming a factor of 2 drop off/increase - the balance drops off with factor 4.
    '''    
    populations = np.zeros(len(ionization_rates)+1)
    #print(ionization_rates)

    populations[0] = 1.0
    for z in range(0,len(ionization_rates)):
        populations[z+1] = populations[z] * ionization_rates[z] / recombination_rates[z]
    
    pp = populations[-1]
    norm = populations.sum()
    current = pp 
    totalarray = [*populations]
    for z in range(len(ionization_rates),atomicNumber):
        current = current * 0.25 
        norm += current 
        totalarray.append(current)
        
    
    totalarray = np.array(totalarray)
    totalarray /= norm
    populations = populations / norm
    #print(populations)
    return populations