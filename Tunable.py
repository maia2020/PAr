"""Code basé sur l'article "Tunable bandgaps in a deployable metamaterial" pour reproduire la matrice de transfert"""


import numpy as np
import matplotlib.pyplot as plt
import sympy as sy

omega=100
rho=7850
L=0.1
E=2.1 * 10**(9)
I=8.3333 * 10**(-13)
k=(( (rho*A) / (E*I) ) * omega**2)**(1/4)
kt=4.2
ligne1 = [  (np.cos(L*k) + np.cosh(L*k))/2  , (np.sin(L*k) + np.sinh(L*k))/2*k , -(np.sin(L*k) - np.sinh(L*k))/2*E*I*(k**3) , (np.cos(L*k) - np.cosh(L*k))/2*E*I*(k**2) ]
ligne2 = [-k*(np.sin(L*k) - np.sinh(L*k))/2  , (np.cos(L*k) + np.cosh(L*k))/2 , (np.cos(L*k) - np.cosh(L*k))/2*E*I*(k**2) , -(np.sin(L*k) + np.sinh(L*k))/2*E*I*k ]
ligne3 = [E*I*(k**3)*(np.sin(L*k) + np.sinh(L*k))/2  , -E*I*(k**2)*(np.cos(L*k) - np.cosh(L*k))/2 , (np.cos(L*k) + np.cosh(L*k))/2 , k*(np.sin(L*k) - np.sinh(L*k))/2 ]
ligne4 = [E*I*(k**2)*(np.cos(L*k) - np.cosh(L*k))/2 , E*I*k*(np.sin(L*k) - np.sinh(L*k))/2 , -(np.sin(L*k) + np.sinh(L*k))/2*k , (np.cos(L*k) + np.cosh(L*k))/2  ]
Transfer=np.array([ligne1,ligne2,ligne3,ligne4]) 


n=0
kv=4
mr= 10
t=1


Tk=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[((n*kv*mr*(omega**2))*np.exp(np.sqrt(-L*omega*t))/((n*kv)-(mr*(omega**2)))),0,0,1]])

T = Transfer@Tk
print(T)

""" print(ligne1)
print(ligne2) """