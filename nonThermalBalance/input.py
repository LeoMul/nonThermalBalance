class input:
    '''
    Class for input, read from json.
    '''
    def __init__(self,
                 #data and elements
                 listOfAtomicNumbers=[],        
                 pathsOfRecombinationData=[],
                 massesOfElements=[],
                 maxIonizationPlus      = 6,
                 
                 #for scaling deposition to the entire ejecta if necessary.
                 massAllEjecta = None,     
                 averageAtomicMass  = 140,
                 
                 # plasma parameters
                 thermalElectronTemperature = None, 
                 imposedElectronDensitySF    = None,
                 imposedElectronDensityRecombination    = None,
                 
                 # transient properties.
                 velocityExpansionC = None,         
                 timeSinceExplosionDays = None, 
                 
                 #deposition:
                 depositionMode = "artis",
                 depfactor      = None,  #override this deposition by factor...
                 depositionOverride = None, #override this deposition to a number
                 
                 #Should I iterate?
                 selfConsistent = False,
                 
                 #photon recycling
                 photonRecycling = False,
                 phi_r           = 0.0, 
                 
                 #Kasen-Barnes paramter
                 velocityMaxForEfficiency   = 0.1,
                 
                 

                 ):        

        #Boring transfer of memory. 
        self.listOfAtomicNumbers        = listOfAtomicNumbers
        self.massesOfElements           = massesOfElements
        self.pathsOfRecombinationData   = pathsOfRecombinationData
        self.thermalElectronTemperature = thermalElectronTemperature
        self.imposedElectronDensitySF   = imposedElectronDensitySF
        self.imposedElectronDensityRecombination = imposedElectronDensityRecombination
        self.velocityExpansionC         = velocityExpansionC
        self.timeSinceExplosionDays     = timeSinceExplosionDays
        self.averageAtomicMass          = averageAtomicMass
        self.maxIonizationPlus              = maxIonizationPlus
        self.selfConsistent             = selfConsistent
        self.depositionOverride         = depositionOverride
        self.velocityMaxForEfficiency   = velocityMaxForEfficiency
        self.depositionMode             = depositionMode.lower() 
        self.depfactor                  = depfactor
        self.photonRecycling            = photonRecycling
        self.phi_r                      = phi_r
        self.massAllEjecta              = massAllEjecta