"""
Eric Figueroa
Ch460 
Professor Topper
Homework#7: Code for Questions 3, 4, 5, and 6. 
"""
import numpy as np # Need numpy for ln, exp, pi, and array functionality
import matplotlib.pyplot as plt # Need matplotlib for plots

k = 1.380658E-23 # Boltzman constant in J/K
Na = 6.0221367E+23 # Avogadro's number
Eh = 4.35974820E-18 # 1 Hartree = 4.35974820E-18 J
a0 = 5.291772109E-11 # Bohr radius in meters
Ang = 1E-10 # 1 Angstrom in meters
# Define paramters for the potential V(x)
E0 = 2E-4           # In Hartrees
C = 0.045           # Unitless
Vh = E0 / (4*C)     # Height of potential barrier 
Xmin = (2*C)**-0.5  # Location of righthand minima

######################## Question 3 ####################################
V = lambda x: E0*(C*x**4 - x**2) # Define the potential
x = np.linspace(-6,6,10000)  # Array of x values 
Vplot = V(x)*Eh*Na/1000     # Potential converted to kJ/mol
Xplot = x*a0/Ang            # X values converted to Angstroms

# Plot V(x) as function of x
plt.figure("V(x)")
plt.plot(Xplot , Vplot,'-r') 
plt.xlabel("x in Angstorms")           
plt.ylabel("Potential V(x) in kJ/mol")
plt.title("Double-well potential as a Function of x")
plt.grid(True, which = 'major', axis = 'both')

######################## Question 4 ####################################
# Redefine the potential such that V at minima is the zero point energy
V = lambda x: (E0*(C*x**4 - x**2) + Vh)*Eh # Converted to J

# Define the distribution as a function of x and temperature T: F(x,T)
F = lambda x,T: np.exp(- V(x) / (k*T)) # T in K and x in atomic units

T0 = 50   # Define value of T that localizes F(x,T) to the minima
T1 = 500  # Define value of T that allows F to spread across both minima 

Fplot = [F(x, T0), F(x, T1)] # List of F trapped and untrapped arrays

plt.figure("Ftrap") # Plot trapped F
plt.plot(Xplot , Fplot[0],'-b') 
plt.xlabel("x in Angstorms")           
plt.ylabel("F(x)")
plt.title("F(x) at 50 K")
plt.grid(True, which = 'major', axis = 'both')
plt.figure("Funtrapped") # Plot untrapped F
plt.plot(Xplot , Fplot[1],'-g')
plt.xlabel("x in Angstorms")           
plt.ylabel("F(x)")
plt.title("F(x) at T = 500 K")
plt.grid(True, which = 'major', axis = 'both')

######################## Question 5/6 ####################################
A = 6       # End of x interval in atomic units
n = 100000   # Number of points in interval 

# Code crude 1D MonteCarlo integrtation for a function f(x,T)
# Input function f(), start of interval, end point, number of points, and T
# Output integral value and standard error of the integral 
def MonteCarlo_Integral(f, a, b, m, T): 
    flist = np.array([])        # Store values of f(x) in array
    fsqrlist = np.array([])     # Store values of f^2 in array
    
    for i in range(0,m): # Generate m random points
        xrandom = a + (b-a)*np.random.rand()       # Random x in interval (a,b)
        frandom = f(xrandom, T)                    # Calculate f at random x
        flist = np.append(flist, frandom)          # Add f to array
        fsqrlist = np.append(fsqrlist, frandom**2)    # Add f^2 to array

    # Calculate the MC integral from the average of f points
    avgf = np.average(flist)
    MCintegral = (b-a)*avgf
    # Calculate the error associated with the MC integral
    avgfsqr = np.average(fsqrlist)
    MCstandarderror = np.sqrt((avgfsqr - avgf**2) / (m - 1)) * (b-a)
    
    return np.array([MCintegral, MCstandarderror])

# Calculate crude MC integrals of F(x,T0) and F(x,T1) in interval (-A,A) 
FMC = [MonteCarlo_Integral(F, 0, A, n, T0), MonteCarlo_Integral(F, 0, A, n, T1)]
Fintegral0 = 2*FMC[0][0] # MC integral is doubled due to symmetry
Fintegral1 = 2*FMC[1][0]  
FMCerror0 = 2*FMC[0][1]  # Error is doubled to match magnitude of Fintegral
FMCerror1 = 2*FMC[1][1]  
print("The integral of F(x, T0) is:", Fintegral0, "±", FMCerror0)
print("The integral of F(x, T1) is:", Fintegral1, "±", FMCerror1, "\n")

# Define G(x) as a function of x and temperature T
G = lambda x,T: V(x)*F(x,T) # T in K and x in atomic units

# Run crude MC algorithm with n points to calculate integral in (-A,A)
GMC = [MonteCarlo_Integral(G, 0, A, n, T0), MonteCarlo_Integral(G, 0, A, n, T1)]
Gintegral0 = 2*GMC[0][0] # MC integral is doubled due to symmetry
Gintegral1 = 2*GMC[1][0]  
GMCerror0 = 2*GMC[0][1]  # Error is doubled to match magnitude of Gintegral
GMCerror1 = 2*GMC[1][1]  
print("The integral of G(x, T0) is:", Gintegral0, "±", GMCerror0)
print("The integral of G(x, T1) is:", Gintegral1, "±", GMCerror1, "\n")

# Compute the average potential energy at each temperature
Vaverage0 = Gintegral0 / Fintegral0
Vaverage1 = Gintegral1 / Fintegral1
print("<V> at", T0,"K is approximately:", Vaverage0, "J")
print("<V> at", T1,"K is approximately:", Vaverage1, "J")

# Compute the average thermodynamic energy in J/mol 
Eaverage0 = (k*T0/2 + Vaverage0) * Na
Eaverage1 = (k*T1/2 + Vaverage1) * Na
print("<E> at", T0,"K is approximately:", Eaverage0, "J / mol")
print("<E> at", T1,"K is approximately:", Eaverage1, "J / mol")

# print("Actual error0 was:",Fintegral0-2.30458) 
# print("Actual error1 was:",Fintegral1-7.5234, "\n") 
