# Ch460 Stat Mech HW#6
""""" Problem 9.3 """""
from numpy import pi
from numpy import exp
N = 6.022e23
kB = 1.3806e-23
h = 6.626e-34
R = 8.314   # Assume 1 mol
P = 101325 # Pascal; Assume at atmospheric pressure
T = [800 + 273, 900 + 273, 1000 + 273, 1100 + 273, 1200 + 273]
Kp_exp = [1.14e-2, 4.74e-2, 0.165, 0.492, 1.23]
V = []  # Assume ideal gas
q_I = []
q_I2 = []
for i in range (0, 5):
    vol = (R*T[i])/P
    V.append(vol)
Kp_calc = []
w = 214 # cm^-1
T_vib = 308 # K
B = 0.0373 # cm^-1
T_rot = 0.0537
k = 1.7e-5 # dynes/cm
Do = 35.6*4184/6.0221367E+23  # J/mol
m_I = 126.90447*(1/N)*(1/1000)
m_I2 = 2*m_I
excited1 = (0.94*(1.60218e-19))/N

for i in range (0, 5):
    I = ((2*pi*m_I*kB*T[i])/(h**2))**(3/2)*V[i]*(2 + 2*exp(-excited1/(kB*T[i])))
    q_I.append(I)
    I2 = ((2*pi*m_I2*kB*T[i])/(h**2))**(3/2)*V[i]*(T[i]/(2*T_rot))*((1 - exp(-T_vib/T[i]))**-1)*exp(Do/(kB*T[i]))
    print(exp(Do/(kB*T[i])))
    q_I2.append(I2) 
print("Partition functions",q_I,q_I2,"\n")
for k in range (0, 4):
    K = (P/R)*((q_I[k]**2)/q_I2[k])
    Kp_calc.append(K)
print(Kp_calc)