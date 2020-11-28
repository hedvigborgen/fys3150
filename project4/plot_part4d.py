# Creating arrays for storing needed values for a random matrix
expEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expEnergySquared = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expMagneticMoment_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
totalEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
flips_random = np.zeros((len(temperatures), MCCs))
