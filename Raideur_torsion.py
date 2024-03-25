import numpy as np
import matplotlib.pyplot as plt



# caractéristiques de la poutre
#eta_w=0.0003*217*1j
eta_w=0.001j
E=(1+eta_w)*69e9
#E=69e9      #Module de Young
rho=2710     #masse volumique
b=0.03       #épaisseur
Long=0.595      #Longuer
h=0.015      #hauteur
I=(b*h**3)/12 #Moment d'inertie
EI=E*I       #rigidity
A=h*b        #aire section
m1=Long*b*h*rho #Masse
####

k_HSLDS=1.1423e+04  #Raideur de l'absorbeur TMD (résultat des différents choix de dimensionnement)





N_cell = 2  #Numero cells
N=int(N_cell/2)
freq = np.logspace(1,3,10000) #frequency
phi=0
frf=[]
frf_tmd=[]
Vr=1
l=Long/N_cell
for f in freq:
    s=2*np.pi*f
    T=[]    
    
            
    kap=(( (rho*A) / (E*I) ) * s**2)**(1/4)
    Tii=np.cos(l*kap)+np.cosh(l*kap)
    T12=(np.sin(l*kap)+np.sinh(l*kap))/kap
    T13=-(np.sin(l*kap)-np.sinh(l*kap))/(EI*kap**3)
    T14=(np.cos(l*kap)-np.cosh(l*kap))/(EI*kap**2)
    T21=-(np.sin(l*kap)-np.sinh(l*kap))*kap
    T24=-(np.sin(l*kap)+np.sinh(l*kap))/(EI*kap)
    T31=(EI*kap**3)*(np.sin(l*kap)+np.sinh(l*kap))
    T32=-(EI*kap**2)*(np.cos(l*kap)-np.cosh(l*kap))
    T34=-T21
    T43=-T12
    T42=(EI*kap)*(np.sin(l*kap)-np.sinh(l*kap))
    T11=Tii*np.cos(phi)
    T33=T11


    T.append(0.5*np.array([[T11,T12,T13,T14,np.sin(phi,0)],[T21,Tii,-T14,T24,0,0],[T31,T32,T33,-T21],[-T32,T42,-T12,Tii]]))
    
    kcoef=1/(1/(k_HSLDS+c2*1j*2*np.pi*f)-1/(4*np.pi**2*f**2*m2))
    
    F=np.array([[1,0,0,0],[0,1,0,0],[-kcoef,0,1,0],[0,0,0,1]])
    
    T_p = T[0] @ T[0]
    
    T_tmd=T[0] @ F @ T[0]   
    
    T_tmd = np.linalg.matrix_power(T_tmd, N)
    T_p = np.linalg.matrix_power(T_p, N)
    
    #T=np.dot(T[2],T[1],T[0])
    Vl = Vr/(T_p[2,2] - (T_p[2,3]*T_p[3,2]/T_p[3,3]))        
    Ml = -T_p[3,2]*Vl/(T_p[3,3])
    Yr = (T_p[0,2]*Vl)  +  (T_p[0,3]*Ml)
    
    Vl_tmd = Vr/(T_tmd[2,2] - (T_tmd[2,3]*T_tmd[3,2]/T_tmd[3,3]))        
    Ml_tmd = -T_tmd[3,2]*Vl_tmd/(T_tmd[3,3])
    Yr_tmd = (T_tmd[0,2]*Vl_tmd)  +  (T_tmd[0,3]*Ml_tmd)
    
    auto_va, auto_vetores=np.linalg.eig(T_tmd)

    
    #Yr= -((T[0,2]*T[2,3])/(T[3,2]-(T[2,2]*(T[3,3]/T[2,3]))))  + ((T[0,3]*T[2,2])/((T[2,3]*T[3,2])-(T[2,2]*T[3,3]))) 
    
    frf.append(abs(Yr))
    frf_tmd.append(abs(Yr_tmd))
    
    mu = 1j*np.log(auto_va)/Long
    mu.sort()
    auto_va_list_2.append(mu[2])
    auto_va_list_3.append(mu[3])
    auto_va_list_im_1.append(mu[2].imag)
    auto_va_list_im_2.append(mu[3].imag)
    

plt.loglog(freq,frf,label='Sans absorbeur')
plt.loglog(freq,frf_tmd,label='Avec absorbeur')
plt.title('Frequency Response Function (FRF) of Euler-Bernoulli Beam')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('FRF Amplitude')
plt.legend()
plt.grid(True)
plt.show()

#print(auto_va_list_2)



plt.plot(freq,auto_va_list_2,".")
plt.plot(freq,auto_va_list_3,".")
#plt.plot(freq,auto_va_list_3)
#plt.title('Auto-valeurs de la matrice de transfert')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Constant de propagation')
plt.grid(True)
plt.show()

plt.plot(freq,auto_va_list_im_1,".")
plt.plot(freq,auto_va_list_im_2,".")
#plt.plot(freq,auto_va_list_3)
#plt.title('Auto-valeurs de la matrice de transfert')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Constant de propagation')
plt.grid(True)
plt.show()