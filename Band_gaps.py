#This script is based on the article Local resonator with high-static-low-dynamic stiffness for lowering band
# gaps of flexural wave in beams

import numpy as np
import matplotlib.pyplot as plt
import sympy as sy

rho = 2600
A = 1.602 * (10 ** (-4))
E = 70 * 10 ** (9)
I = 5.968 * (10 ** (-9))
L = 0.125
eta = 1
kv = 1.65 * (10 ** (5))
mr = 0.0437
solutions=[]
omegas_plot=[]
solution=0.5

for omega in range(1,2000):
    omegas_plot.append(omega/(6.28))
    def beta(omega):
        return (( (rho*A) / (E*I) ) * omega**2)**(1/4)

    H1 = [sy.cos(L*beta(omega)), sy.sin(L*beta(omega)), sy.cosh(L*beta(omega)), sy.sinh(L*beta(omega))]
    H2 = [-beta(omega)*sy.sin(L*beta(omega)), beta(omega)*sy.cos(L*beta(omega)), beta(omega)*sy.sinh(L*beta(omega)), beta(omega)*sy.cosh(L*beta(omega))]
    H3 = [-(beta(omega)**2)*sy.cos(L*beta(omega)), -(beta(omega)**2)*sy.sin(L*beta(omega)), (beta(omega)**2)*sy.cosh(L*beta(omega)), (beta(omega)**2)*sy.sinh(L*beta(omega))]
    H4 = [(beta(omega)**3)*sy.sin(L*beta(omega)), -(beta(omega)**3)*sy.cos(L*beta(omega)), (beta(omega)**3)*sy.sinh(L*beta(omega)), (beta(omega)**3)*sy.cosh(L*beta(omega))]
    H = np.array([H1, H2, H3, H4])
    H = sy.Matrix(H)
    #print(H)



    def FG(omega):
        return (1/E*I) * (eta * kv * mr * omega*omega) / ( (eta * kv) - (mr * omega*omega))

    G1 = [1, 0, 1, 0]
    G2 = [0, beta(omega), 0, beta(omega)]
    G3 = [-beta(omega)**2, 0, beta(omega)**2, 0]
    G4 = [-FG(omega), -beta(omega)**3, -FG(omega), beta(omega)**3]
    G = np.array([G1, G2, G3, G4])
    G = sy.Matrix(G)
    #print(G)

    G_inv = sy.Matrix(G).inv()

    Id = np.identity(4)
    q = sy.symbols('q')

    expression = sy.det(G_inv * H - sy.exp(1j * q * L) * Id)
    equation = sy.Eq(expression, 0)

    #solutions = sy.solve(equation, q)
    try:
    # Attempt to find a solution
        solution = sy.nsolve(equation, q, solution)
        solutions.append(solution)
    except:
    # Handle the case where no solution is found
        solutions.append(0)

for i in range(len(solutions)):
    solutions[i]=solutions[i]*L/3.14
plt.plot(omegas_plot,[sy.re(solution) for solution in solutions])
