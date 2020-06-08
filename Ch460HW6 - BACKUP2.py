"""
ERic Figueroa
Ch460 
Professor Topper
Homework#6: Problems 8.5,7 and 13. 
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
print("The entropy of carbon dioxide at 25°C and 1 atm is", S_CO2, "J / K mol \n")

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
    
    # Identify which entropy terms see suspicious to get the right calculation
    # [Na*k*Strans,Na*k*Srot,Na*k*Svib,Na*k*Selec]
    return S 
  
# Methane is a nonlinear molecule with 12-fold symmetry. 
CH4θr = [7.54,7.54,7.54]
CH4θv = [4170,2180,2180,4320,4320,4320, 1870, 1870, 1870]

S_CH4 = S_NonLinear(101325, 298.15, 16.04, CH4θr, CH4θv,12)
print("The entropy of methane at 25°C and 1 atm is", S_CH4, "J / K mol \n")


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
print("The constant volume heat capacity for ammonia at 300 K is ", Cv_NH3, "J / K mol \n")

