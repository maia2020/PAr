'''
Simulation d'un Tuned Mass Damper. Objectifs: Diagramme de Bode en gain et la reponse temporelle
Le seisme qui apparait est une réponse sinusoidale qui represente la fréquence de résonance de la poutre avant de mettre l'absorbeur
basé sur le code de: Gauthier Legrand et Francis Pagaud. 18/06/2020.
'''
import numpy as np # Outils numeriques
import scipy as sp # Outils scientifiques
import scipy.integrate # pour l'integration
import matplotlib.pyplot as plt # Outils graphiques
from matplotlib.widgets import Slider
import matplotlib

#matplotlib.rc('xtick', labelsize=24)
#matplotlib.rc('ytick', labelsize=24)
#matplotlib.rcParams.update({'font.size': 22})


## Parametres et conditions initiales

omega_0 = 1# omega_0^2 = k/m
omega_1 = 1 # omega_0^2 = k_1/m_1
a = 0.01 # a = m_1/m
A0 = 1
eta = 1e-2
j = complex(0,1)
omega = np.logspace(-0.3,0.3,1000) #

def H_1(omega, omega_1, eta) :
    return(omega**2/(-omega**2 + 2j*eta*omega_1*omega + omega_1**2))

def H_2(omega, omega_1, omega_0, a, eta) :
    return(omega**2/((1 + a + a*H_1(omega, omega_1, eta))*omega**2 - omega_0**2))



tfin = 500 # instant final de la simulation
dt = 1e-1 # pas
t = np.linspace(0,tfin,int(tfin/dt))


x_0 = 0
x_10 = 0
dxdt_0 = 0
dx1dt_0 = 0

vec0 = [x_0, x_10, dxdt_0, dx1dt_0]

##Fonctions

#y = [x, x_1, dxdt, dx_1dt]
def systeme(y,t) :
    global omega_0, omega_1, a, eta

    (x, x_1, dxdt, dx_1dt) = y

    d2xdt2 = - omega_0**2*x + 2*eta*omega_1*a*dx_1dt + omega_1**2*a*x_1 + A0*np.sin(omega_0*t)
    d2x_1dt2 =  omega_0**2*x - (1+a)*2*eta*omega_1*dx_1dt - (1+a)*omega_1**2*x_1 - A0*np.sin(omega_0*t)


    res = [dxdt, dx_1dt, d2xdt2, d2x_1dt2]
    return res


def update(val):
    global omega_0, omega_1, a, eta

    eta = 10**(s_eta.val)
    omega_1 = 10**(s_omega1.val)
    a = 10**(s_a.val)


    sol = sp.integrate.odeint(systeme,vec0, t)

    x = sol[:, 0]
    x_1 = sol[:, 1]


    l1.set_ydata(np.abs(H_2(omega, omega_1, omega_0, a, eta)))
    l2.set_ydata(np.abs(H_1(omega, omega_1, eta)))

    lj.set_ydata(x)
    ll.set_ydata(x_1)

    fig.canvas.draw_idle()


## Resolution
 # Determination de lâ€™entree de lâ€™etage dâ€™amplification
sol = sp.integrate.odeint(systeme,vec0,t)

x = sol[:, 0]
x_1 = sol[:, 1]


## Traces des courbes

fig = plt.figure('TMD')
plt.clf()

ax1 = plt.axes([0.07, 0.23, 0.4, 0.7])
ax2 = plt.axes([0.57, 0.23, 0.4, 0.7])

#slider
ax_eta = plt.axes([0.15, 0.1, 0.75, 0.03])
ax_omega1 = plt.axes([0.15, 0.06, 0.75, 0.03])
ax_a = plt.axes([0.15, 0.02, 0.75, 0.03])

s_eta = Slider(ax_eta, '$\log \ \eta$', -3, 0, valinit=np.log10(eta), valfmt='%0.2f')
s_omega1 = Slider(ax_omega1, '$ \log \ \omega_1 / \omega_0$ ', -1, 4, valinit=np.log10(omega_1), valfmt='%0.1f')
s_a = Slider(ax_a, '$\log \ a$', -4, 1, valinit=np.log10(a), valfmt='%0.2f')

s_eta.on_changed(update)
s_omega1.on_changed(update)
s_a.on_changed(update)


#Diagramme de Bode
l1, = ax1.loglog(omega/omega_0, np.abs(H_2(omega, omega_1, omega_0, a, eta)),lw=5, label='$H_2$', color='steelblue') # fonction de transfert de la tour
l2, = ax1.loglog(omega/omega_0, np.abs(H_1(omega, omega_1, eta)),lw=5, label='$H_1$', color='darkorange') # fonction de transfert du TMD


#Reponse dynamique
ll, = ax2.plot(t, x_1, label='$x_1$', lw=3,color= 'darkorange')
lj, = ax2.plot(t, x, label='$x$', lw=3, color='steelblue')

ax1.set_xlabel('$\omega / \omega_0$', fontsize = 28)
ax1.set_ylabel('$|H|$', fontsize = 28)
ax1.set_ylim(10**(-1), 10**2)
ax1.legend(framealpha = 0.5, loc=2)
ax1.grid()
#ax2.set_ylim(-0.05, 1.5)
ax2.set_xlabel("$t$ (s)")
ax2.set_ylabel("RÃ©ponse dynamique")
ax2.legend(framealpha = 0.5, loc=2)


mng = plt.get_current_fig_manager()     #Plein ecran

plt.show()