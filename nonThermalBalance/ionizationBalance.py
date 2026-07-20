import numpy as np 
def ionizationBalance(ionization_rates, recombination_rates,atomicNumber):
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
    
    return populations