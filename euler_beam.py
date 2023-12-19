import numpy as np
import matplotlib.pyplot as plt


# caractéristiques de la poutre

E=69e9      #Module de Young
#eta_w=0.0003*217*1j
E=(1+0.01j)*69e9      #Module de Young
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



#lfreq=[34.766 217.01 596.4]
freq=34.766         

ksi1=0.01        
w1=2*np.pi*freq 
k1=(w1**2)*m1 
c1=2*ksi1*w1*m1

mu=0.05
ksi2op=np.sqrt(3*mu/8/(1+mu))+0.1616*ksi1/(1+mu)
qop=1/(1+mu)*(1-1.5906*ksi1*np.sqrt(mu/(1+mu)))

w2=qop*w1
m2=mu*m1
c2=2*ksi2op*w2*m2
k_HSLDS=(w2**2)*m2



N_cell = 2  #Numero cells
freq = np.logspace(1,3,10000) #frequency
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


    T.append(0.5*np.array([[Tii,T12,T13,T14],[T21,Tii,-T14,T24],[T31,T32,Tii,-T21],[-T32,T42,-T12,Tii]]))
    
    kcoef=1/(1/(k_HSLDS+c2*1j*2*np.pi*f)-1/(4*np.pi**2*f**2*m2))
    
    F=np.array([[1,0,0,0],[0,1,0,0],[-kcoef,0,1,0],[0,0,0,1]])
    
    T_p = T[0] @ T[0]
    
    T_tmd=T[0] @ F @ T[0]   
    #T = np.linalg.matrix_power(T, N_cell)
    
    #T=np.dot(T[2],T[1],T[0])
    Vl = Vr/(T_p[2,2] - (T_p[2,3]*T_p[3,2]/T_p[3,3]))        
    Ml = -T_p[3,2]*Vl/(T_p[3,3])
    Yr = (T_p[0,2]*Vl)  +  (T_p[0,3]*Ml)
    
    Vl_tmd = Vr/(T_tmd[2,2] - (T_tmd[2,3]*T_tmd[3,2]/T_tmd[3,3]))        
    Ml_tmd = -T_tmd[3,2]*Vl_tmd/(T_tmd[3,3])
    Yr_tmd = (T_tmd[0,2]*Vl_tmd)  +  (T_tmd[0,3]*Ml_tmd)
    
    
    
    
    #Yr= -((T[0,2]*T[2,3])/(T[3,2]-(T[2,2]*(T[3,3]/T[2,3]))))  + ((T[0,3]*T[2,2])/((T[2,3]*T[3,2])-(T[2,2]*T[3,3]))) 
    
    frf.append(abs(Yr))
    frf_tmd.append(abs(Yr_tmd))
    

plt.loglog(freq,frf)
plt.loglog(freq,frf_tmd)
plt.title('Frequency Response Function (FRF) of Euler-Bernoulli Beam')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('FRF Amplitude')
plt.grid(True)
plt.show()