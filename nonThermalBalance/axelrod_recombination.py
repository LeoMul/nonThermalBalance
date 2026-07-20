import numpy as np 

'''
Recombination rate coefficients, for when the user does not supply recombination rates.
'''

def RRAxelrodGround(temp_kelvin,i):
    #i = ionization
    
    temp_1e4K = temp_kelvin / 1e4 
    
    return 1E-13 * i*i * np.power(temp_1e4K,-0.5 )

def RRAxelrodTotal(temp_kelvin,i):
    #i = ionization
    
    temp_1e4K = temp_kelvin / 1e4 
    
    return 3E-13 * i*i * np.power(temp_1e4K,-0.75 )

