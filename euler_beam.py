import numpy as np
import matplotlib.pyplot as plt


# caractéristiques de la poutre

E=69e9      #Module de Young
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
freq = np.logspace(1,3,10000) #frequency
frf=[]
Vr=1

for f in freq:
    s=2*np.pi*f
    coef=[0.5, 0.499 ,0.001]
    l1=Long*coef[0]
    l2=Long*coef[1]
    l3=Long*coef[2]
    l=[l1, l2, l3]
    T=[]    
    for cell in range(3):
            
        kap=(( (rho*A) / (E*I) ) * s**2)**(1/4)
        Tii=np.cos(l[cell]*kap)+np.cosh(l[cell]*kap)
        T12=(np.sin(l[cell]*kap)+np.sinh(l[cell]*kap))/kap
        T13=-(np.sin(l[cell]*kap)-np.sinh(l[cell]*kap))/(EI*kap**3)
        T14=(np.cos(l[cell]*kap)-np.cosh(l[cell]*kap))/(EI*kap**2)
        T21=-(np.sin(l[cell]*kap)-np.sinh(l[cell]*kap))*kap
        T24=-(np.sin(l[cell]*kap)+np.sinh(l[cell]*kap))/(EI*kap)
        T31=(EI*kap**3)*(np.sin(l[cell]*kap)+np.sinh(l[cell]*kap))
        T32=-(EI*kap**2)*(np.cos(l[cell]*kap)-np.cosh(l[cell]*kap))
        T34=-T21
        T43=-T12
        T42=(EI*kap)*(np.sin(l[cell]*kap)-np.sinh(l[cell]*kap))

        T.append(0.5*np.array([[Tii,T12,T13,T14],[T21,Tii,-T14,T24],[T31,T32,Tii,-T21],[-T32,T42,-T12,Tii]]))
        
        

    
    T=T[2]@T[1]@T[0]
    T=T/1000
    print(T)

    #T=np.dot(T[2],T[1],T[0])
    Vl = Vr/(T[2,2] - (T[2,3]*T[3,2]/T[3,3]))        
    Ml = -T[3,2]*Vl/(T[3,3])
    Yr = (T[0,2]*Vl)  +  (T[0,3]*Ml)
    
    #Yr= -((T[0,2]*T[2,3])/(T[3,2]-(T[2,2]*(T[3,3]/T[2,3]))))  + ((T[0,3]*T[2,2])/((T[2,3]*T[3,2])-(T[2,2]*T[3,3]))) 
    
    frf.append(abs(Yr))
    

plt.loglog(freq, frf)
plt.title('Frequency Response Function (FRF) of Euler-Bernoulli Beam')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('FRF Amplitude')
plt.grid(True)
plt.show()

        
    