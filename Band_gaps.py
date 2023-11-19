#This script is based on the article Local resonator with high-static-low-dynamic stiffness for lowering band
# gaps of flexural wave in beams

import numpy as np
import matplotlib.pyplot as plt
import sympy as sy

rho=2600
A=1.602*(10**(-4))
E=70*10**(9)
I=5.968*(10**(-9))
omega = np.linspace(1, 1000, 2000)
L=0.125

def beta(omega):
    return (( (rho*A) / (E*I) ) * omega**2)**(1/4)


H1 = [  np.cos(L*beta(omega))  , np.sin(L*beta(omega)), np.cosh(L*beta(omega)) , np.sinh(L*beta(omega)) ]
H2 = [  -beta(omega)*np.sin(L*beta(omega))  , beta(omega)*np.cos(L*beta(omega)), beta(omega)*np.sinh(L*beta(omega)) , beta(omega)*np.cosh(L*beta(omega)) ]
H3 = [  -(beta(omega)**2)*np.cos(L*beta(omega))  , -(beta(omega)**2)*np.sin(L*beta(omega)), (beta(omega)**2)*np.cosh(L*beta(omega)) , (beta(omega)**2)*np.sinh(L*beta(omega)) ]
H4 = [  (beta(omega)**3)*np.sin(L*beta(omega))  , -(beta(omega)**3)*np.cos(L*beta(omega)), (beta(omega)**3)*np.sinh(L*beta(omega)) , (beta(omega)**3)*np.cosh(L*beta(omega)) ]
H=np.array([H1,H2,H3,H4]) 

eta=1
kv=1.65*(10**(5))
mr=0.0437

def FG(omega):
    return (1/E*I) * (eta * kv * mr * omega*omega) / ( (eta * kv) - (mr * omega*omega))

G1=[1,0,1,0]
G2=[0,beta(omega),0,beta(omega)]
G3=[-beta(omega)**2,0,beta(omega)**2,0]
G4=[-FG(omega),-beta(omega)**3,-FG(omega),beta(omega)**3]
G=np.array([G1,G2,G3,G4]) 

G_inv = np.linalg.inv(G)

Id=np.identity(4)
q = sy.symbols('q')

expression = sy.det(G_inv * H - sy.exp( 1j * q * L) * Id)
equation = sy.Eq(expression, 0)

solutions = sy.solve(equation, q)

print("Solutions for q:")
for solution in solutions:
    print(solution)