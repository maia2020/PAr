"""
Cette première partie a comme objectif de restituer les resultats de l'article "Vibration isolation characteristics of a nonlinear isolator
using Euler buckled beam as negative stiffness corrector: A theoretical and experimental study" de Xiuchang Huang et al.
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sy


# Paramètres du système
q0tilda = 0.03  # non dimensional vertical displacement
theta = 25  # Lateral deflection angle
gamma = np.cos(theta * np.pi / 180)

# On veut tracer la courbe de la figure 3 en resolvant l'equation (3) de l'article

a = np.sqrt(((np.pi**2) * q0tilda**2) - 4 * gamma + 4)
b = np.pi * q0tilda

k1 = (a - b / (a * gamma)) * (((b**2) / 2) - 2 * gamma + 6)
k3 = (
    a
    - (b / (a * gamma**2))
    + (a - (b / (2 * a * gamma**3)) + (b / ((gamma**2) * a**3)))
    * ((b**2) / 2 - 2 * gamma + 6)
)


# Alors on a l'equation suivante:
def Ftilta(u):
    return -k1 * u + k3 * u**3


# On crée un vecteur de valeurs de utilda
utilda = np.linspace(-0.4, 0.4, 500)

# On crée un vecteur de valeurs de Ftilta
Ftilta = Ftilta(utilda)

# On trace la courbe
plt.plot(utilda, Ftilta)
plt.xlabel("utilda")
plt.ylabel("Ftilta")
plt.title("Figure 3")
plt.show()

# Alors on va différencier Ftilta pour trouver la raideur
# On va utiliser python pour faire la derivation symbolique

utilda = sy.Symbol("utilda")
gamma = sy.Symbol("gamma")
q0 = sy.Symbol("q0")

Ftilda = (
    (
        1
        - sy.pi
        * q0
        * (sy.pi**2 * q0**2 + 4 * (1 - sy.sqrt(utilda**2 + gamma**2)))
        ** (-1 / 2)
    )
    * (2 * sy.sqrt(utilda**2 + gamma**2) - (12 + sy.pi**2 * q0**2) / 2)
    * (utilda / sy.sqrt(utilda**2 + gamma**2))
)

# On fait la derivation
raideurNomDim = sy.diff(Ftilda, utilda)

# Alors on plot la raideur en fonction de utilda

raideurNomDim = sy.lambdify((utilda, gamma, q0), raideurNomDim)

q0tilda = 0.03  # non dimensional vertical displacement
theta = 25  # Lateral deflection angle
gamma = np.cos(theta * np.pi / 180)
# On crée un vecteur de valeurs de utilda
utilda = np.linspace(-0.4, 0.4, 500)

# On crée un vecteur de valeurs de raideur
raideur = raideurNomDim(utilda, gamma, q0tilda)

# On trace la courbe
plt.plot(utilda, raideur)
plt.xlabel("utilda")
plt.ylabel("raideur")
plt.title("Raideur en fonction de utilda")
plt.show()
