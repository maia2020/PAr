"""Code basé sur l'article "Tunable bandgaps in a deployable metamaterial" pour reproduire la matrice de transfert"""


import numpy as np
import matplotlib.pyplot as plt
import sympy as sy

L=1
k=1
E=1
I=1
ligne1 = [  (np.cos(L*k) + np.cosh(L*k))/2  , (np.sin(L*k) + np.sinh(L*k))/2*k , -(np.sin(L*k) - np.sinh(L*k))/2*E*I*(k**3) , (np.cos(L*k) - np.cosh(L*k))/2*E*I*(k**2) ]
ligne2 = [-k*(np.sin(L*k) - np.sinh(L*k))/2  , (np.cos(L*k) + np.cosh(L*k))/2 , (np.cos(L*k) - np.cosh(L*k))/2*E*I*(k**2) , -(np.sin(L*k) + np.sinh(L*k))/2*E*I*k ]
ligne3 = [E*I*(k**3)*(np.sin(L*k) + np.sinh(L*k))/2  , -E*I*(k**2)*(np.cos(L*k) - np.cosh(L*k))/2 , (np.cos(L*k) + np.cosh(L*k))/2 , k*(np.sin(L*k) - np.sinh(L*k))/2 ]
ligne4 = [E*I*(k**2)*(np.cos(L*k) - np.cosh(L*k))/2 , E*I*k*(np.sin(L*k) - np.sinh(L*k))/2 , -(np.sin(L*k) + np.sinh(L*k))/2*k , (np.cos(L*k) + np.cosh(L*k))/2  ]
Transfer=np.array([ligne1,ligne2,ligne3,ligne4]) 

cell_gauche=np.array([1,1,1,1]).T


cell_droit = Transfer@cell_gauche
print(cell_droit)
""" print(ligne1)
print(ligne2) """