"""
Eric Figueroa
Ch460 
Professor Topper
Homework#6: Code for Problems 8.5, 8.7, 8.13, 9.3, 9.8, and 9.9
"""

# Need numpy for column vectors
import numpy as np

# Need pyplot from the matplotlib module to plot
import matplotlib.pyplot as plt

# Define important constants: 
k = 1.380658E-23 # Boltzman constant in J/K
Na = 6.0221367E+23 # Avogadro's number
h = 6.6260755E-34 # Planck's constant in J*s
amu_conversion = 1.6605402E-27 # Convert amu to kg
Pa2atm = 101325 # Pascals in 1 atm

# Problem 8.5: 

# Molar Entropy for an ideal linear polyatomic molecule based off Equation 8-27.
# Takes in inputs of Pressure (P) in Pa and Molecular mass (M) in kg/mol
# Also takes inputs of Temperature (T), Rotational Temperature (θr),
# and Vibrational Temperature (θv) in Kelvin.  θv must be an array.  
# Also inputs Rotational Symmetry number as sym.
def S_Linear(P, T, M, θr, θv, sym): 
    
    We1 = 1 # Ground electronic state degeneracy is usually 1: Change if not.  
    n = len(θv) # Number of vibrational frequencies
    M = M*amu_conversion # Convert molecular mass (amu) to mass in kilograms 
    
    # Calculate the volume per number of molecules: 
    V = k*T/P

    # Solve for the first term (Translational contibution to S)
    Strans = np.log((((2*np.pi*M*k*T)/(h**2))**1.5)*(V*np.e**2.5))
    
    # Solve for the second term (Rotational contribution to S)
    Srot = np.log((T*np.e)/(sym*θr))
    
    # Solve for the third term (Vibrational contribution to S)
    Svibterms = [] # Array for vibrational terms
    
    for i in range(n):
        # Calculate term in the vibrational entropy summation
        term = (θv[i]/T)/(np.exp(θv[i]/T)-1) - np.log(1-np.exp(-θv[i]/T))
        
        Svibterms.append(term) # Add term to array

    Svib = np.sum(Svibterms) # Compute the summation of vibrational terms
    
    # Solve for the fourth entropy term (Electronic contribution to S)
    Selec = np.log(We1) # Will be 0 for We1 = 1
                   
    # Solve for the total molar entropy and return output
    S = Na*k*(Strans + Srot + Svib + Selec)
    
    # Identify which entropy terms see suspicious to get the right calculation
    # [Na*k*Strans,Na*k*Srot,Na*k*Svib,Na*k*Selec]
    return S 

# Carbon dioxide is a linear molecule with a center of Symmetry, sym = 2
S_CO2 = S_Linear(101325, 298.15, 44.01, 0.561, [3360, 954, 954, 1890],2)
print("The entropy of carbon dioxide at 25°C and 1 atm is", S_CO2, "J/K mol \n")

# Problem 8.7: 

# Molar Entropy for an ideal nonlinear polyatomic molecule based off Equation 8-33.
# Takes in inputs of Pressure (P) in Pa and Molecular mass (M) in kg/mol
# Also takes inputs of Temperature (T), Rotational Temperature (θr),
# and Vibrational Temperature (θv) in Kelvin. θr and θv must be arrays.  
# Also inputs Rotational Symmetry number as sym.
def S_NonLinear(P, T, M, θr, θv, sym): 
    
    We1 = 1 # Ground electronic state degeneracy is usually 1: Change if not.  
    n = len(θv) # Number of vibrational frequencies
    M = M*amu_conversion # Convert molecular mass (amu) to mass in kilograms 
    
    # Calculate the volume per number of molecules: 
    V = k*T/P

    # Solve for the first term (Translational contibution to S)
    Strans = np.log((((2*np.pi*M*k*T)/(h**2))**1.5)*(V*np.e**2.5))
    
    # Solve for the second term (Rotational contribution to S)
    Srot = np.log(((np.pi**0.5)*(np.e**1.5)/sym)*((T**3)/(θr[0]*θr[1]*θr[2]))**0.5)
    
    # Solve for the third term (Vibrational contribution to S)
    Svibterms = [] # Array for vibrational terms
    
    for i in range(n):
        # Calculate each term in the vibrational entropy summation
        term = (θv[i]/T) / (np.exp(θv[i]/T)-1) - np.log(1-np.exp(-θv[i]/T))
        
        Svibterms.append(term) # Add term to array

    Svib = np.sum(Svibterms) # Compute the summation of vibrational terms
    
    # Solve for the fourth entropy term (Electronic contribution to S)
    Selec = np.log(We1) # Will be 0 for We1 = 1
                   
    # Solve for the total molar entropy and return output
    S = Na*k*(Strans + Srot + Svib + Selec)
    return S 
  
# Methane is a nonlinear molecule with 12-fold symmetry. 
CH4θr = [7.54,7.54,7.54]
CH4θv = [4170,2180,2180,4320,4320,4320, 1870, 1870, 1870]

S_CH4 = S_NonLinear(101325, 298.15, 16.04, CH4θr, CH4θv,12)
print("The entropy of methane at 25°C and 1 atm is", S_CH4, "J/Kmol \n")

# Problem 8.13: 

# Cv for an ideal nonlinear polyatomic molecule based off Equation 8-32.
# Takes inputs of Temperature (T) and Vibrational Temperature (θv) in Kelvin. 
# θv must be an array. 
def Cv_NonLinear(T, θv): 

    n = len(θv) # Number of vibrational frequencies

    # Calculate the contribution of vibrational motion to the heat capacity Cv
    Cvvibterms = [] # Array for vibrational terms
    
    for i in range(n):
        # Calculate each term in the summation
        term = ((θv[i]/T)**2) * np.exp(θv[i]/T) / ((np.exp(θv[i]/T)-1)**2)
        Cvvibterms.append(term) # Add term to array

    Cvvib = np.sum(Cvvibterms) # Compute the summation of vibrational terms

    # Solve for the total molar heat capacity: 
    # The two 3/2 are the translational and rotational contributions to Cv
    Cv = Na*k*(3/2 + 3/2 + Cvvib)
    return Cv

Cv_NH3 = Cv_NonLinear(300, [4800,1360,4880,4880,2330,2330])
print("The constant volume heat capacity for ammonia at 300 K is", Cv_NH3, "J/Kmol \n")

# Question 9.3:

# Calculate the molecular partition function for an ideal monatomic gas per meter cubed
# Takes inputs of Temperature (T), Molecular mass (M) in amu, and electronic degeneracy (We1)
def q_Monatomic_per_m3(T, M, We1):
    
    M = M*amu_conversion # Convert molecular mass (amu) to mass in kilograms 
    
    # Calculate translational partition function
    qtrans = (((2*np.pi*M*k*T)/(h**2))**1.5)
    
    # Calculate electronic partition function: modify if excited states are accessible
    qelec = We1
    
    # Calculate the partition function per meter cubed
    q = qtrans * qelec
    return q

# Calculate the molecular partition function for an ideal diatomic gas per meter cubed
# Input Temperature (T) in K, Molecular mass (M) in amu, Vibrational Temperature (θv) in K,
# Rotational temperature (θr) in K, and energy of dissociation (Do) in J/mol
def q_Diatomic_per_m3(T, M, θr, θv, Do):
    
    We1 = 1 # Ground electronic state degeneracy is usually 1: Change if not.  
    M = M*amu_conversion # Convert molecular mass (amu) to mass in kilograms 
    Do = Do / Na # Convert Do from J / mol to J
    
    # Calculate translational partition function
    qtrans = (((2*np.pi*M*k*T)/(h**2))**1.5)
    
    # Calculate rotational partition function
    qrot = T/(2*θr)
    
    # Calculate vibrational partition function
    qvib = (1 - np.exp(-θv/T))**-1
    
    # Calculate electronic partition function: modify if excited states are accessible
    qelec = We1 * np.exp(Do / (k*T))
    #print(qelec)
    
    # Calculate the partition function per meter cubed
    q = qtrans * qrot * qvib * qelec
    return q

# Calculate the equilibrium constant Kp(T) using the stoichemetry and the
# partition functions of the products and reactants per meter cubed
# Input Temperature (T) in K. If no species B or D, input 0 for v and 1 for q
def  Kp_Const(T, qA, qB, qC, qD, vA, vB, vC, vD): 
    
    # Calculate the contributions made by concentrations of each species
    Kc = (qC**vC) * (qD**vD) * (qA**-vA) * (qB**-vB)
    
    # Calculate Kp
    Kp = ((k*T)**(vC + vD - vA - vB)) * Kc
    return Kp

# Calculate the equilibrium constant Kp from T = 800 K to T = 1200 K
T = [800,900,1000,1100,1200,872, 973, 1073, 1173, 1274]

for i in range(0, len(T)):
    # Calculate the partition function per meter cubed for monatomic iodine at T
    qI = q_Monatomic_per_m3(T[i], 126.90447, 4)
    
    # Calculate the partition function per meter cubed for diatomic iodine at T
    qI2 = q_Diatomic_per_m3(T[i], 253.8089, 0.0537, 308, 35.6*4184)
    
    #print([qI,qI2])
    # Calculate Kp at T. Input partitions functions and stoichemetry 
    Kp_iodine = Kp_Const(T[i], qI2, 1, qI, 1, 1, 0, 2, 0)
    print("At",T[i], ", Kp = ", Kp_iodine / Pa2atm)
    
print("For the reaction I2 <-> 2I in 1/atm \n")

# DEBUG WITH Na
# qNa = q_Monatomic_per_m3(1000, 22.989769, 2)
# qNa2 = q_Diatomic_per_m3(1000, 22.989769*2, 0.221, 229, 17.3*4180)
# print("Na DEBUG:")
# print("Reactant (Monatomic):", qNa)
# print("Product (Diatomic):", qNa2)
# print("Kp:", Kp_Const(1000, qNa, 1, qNa2, 1, 2, 0, 1, 0)*Pa2atm)
# print(" \n")

# Question 9.9

# Calculate the equilibrium constant Kp from T = 1500 K to T = 2500 K
T = [1500, 2000, 2500]

for i in range(0, len(T)):
    # Calculate the partition function per meter cubed for diatomic nitrogen
    qN2 = q_Diatomic_per_m3(T[i], 28.006148, 2.88, 3374, 225.1*4184)
    
    # Calculate the partition function per meter cubed for diatomic oxygen 
    qO2 = q_Diatomic_per_m3(T[i], 31.9988, 2.07, 2256, 118*4184)
    
    # Calculate the partition function per meter cubed for nitric oxide
    qNO = q_Diatomic_per_m3(T[i], 30.002474, 2.45, 2719, 150*4184)
    
    #print([qN2,qO2,qNO])
    
    # Calculate Kp at T. Input partitions functions and stoichemetry 
    Kp_nitric = Kp_Const(T[i], qN2, qO2, qNO, 1, 0.5, 0.5, 1, 0)
    print("At",T[i], ", Kp = ", Kp_nitric)
    
print("For the reaction N2 + O2 <-> NO \n")

# Question 9.8

# Calculate the equilibrium constant Kc(T) using the stoichemetry and the
# partition functions of the products and reactants per meter cubed
# If no species B or D, input 0 for v and 1 for q(T)
def  Kc_Const(qA, qB, qC, qD, vA, vB, vC, vD): 
    
    Kc = (qC**vC) * (qD**vD) * (qA**-vA) * (qB**-vB)
    return Kc

# Calculate ΔH of a chemical reaction within a short interval [T1, T2]
def ΔH_Reaction(T1, T2, Kc1, Kc2):
    
    ΔH = -Na*k * np.log(Kc2/Kc1) * (T2**-1 - T1**-1)**-1
    return ΔH

T = [300-10**-5,300+10**-5] # Create short temperature range, [T0, T1]
Kc_HI = [] # List of equilibrium constants at T0 and T1

for i in range(0, len(T)):
    # Calculate partition function for I2 at T
    qI2= q_Diatomic_per_m3(T[i], 253.8089, 0.0537, 308, 35.6*4184)
    
    # Calculate partition function for H2 at T
    qH2 = q_Diatomic_per_m3(T[i], 2.01568, 85.3, 6215, 103.2*4184)
    
    # Calculate partition function for HI at T
    qHI = q_Diatomic_per_m3(T[i], 127.91229, 9.06, 3266, 70.5*4184)
    
    # Calculate the equilibrium constants at temperature T[i]: put in list
    Kc_HI.append(Kc_Const(qI2, qH2, qHI, 1, 1, 1, 2, 0))

# Calculate the enthalpy of reaction for the reaction H2 + I2 <-> 2HI
print("At",np.average(T), ", ΔH = ", ΔH_Reaction(T[0], T[1], Kc_HI[0], Kc_HI[1])/4.184)
print("For the reaction H2 + I2 <-> 2HI in cal/mol")


T = np.linspace(300-10**-5,300+10**-5,2) # Create short temperature range about 300 K
Kc_list = Kc_Const(q_Diatomic_per_m3(T, 253.8089, 0.0537, 308, 35.6*4184), # Input q for I2
                    q_Diatomic_per_m3(T, 2.01568, 85.3, 6215, 103.2*4184), # Input q for H2
                    q_Diatomic_per_m3(T, 127.91229, 9.06, 3266, 70.5*4184), # Input q for HI
                    1, 1, 1, 2, 0) # Input stoichemetric coefficients

# Calculate enethalpy of reaction. 
print("At",np.average(T),", ΔH = ",ΔH_Reaction(T[0],T[1],Kc_list[0],Kc_list[1])/4.184)
print("For the reaction H2 + I2 <-> 2HI in cal/mol")

# lnKc = np.log(Kc_list)
# plt.plot(1000/T,lnKc,'r-')

