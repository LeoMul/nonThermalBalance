class input:
    '''
    Class for input, read from json.
    '''
    def __init__(self,
                 listOfAtomicNumbers=[],
                 pathsOfRecombinationData=[],
                 massesOfElements=[],
                 thermalElectronTemperature = None,
                 imposedElectronDensitySF    = None,
                 imposedElectronDensityRecombination    = None,
                 velocityExpansionC = None,
                 timeSinceExplosionDays = None, 
                 selfConsistent = False,
                 averageAtomicMass  = 140,
                 maxIonizationPlus      = 6,
                 depositionOverride = None,
                 velocityMaxForEfficiency   = 0.1,
                 depositionMode = "KasenBarnes",
                 depfactor      = None, 
                 photonRecycling = False,
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