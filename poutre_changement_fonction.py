import numpy as np
import matplotlib.pyplot as plt



# caractéristiques de la poutre
#eta_w=0.0003*217*1j
eta_w=0.001j
E=(1+eta_w)*69e9
#E=69e9      #Module de Young
rho=2710     #masse volumique
b=0.03       #épaisseur
b2=b*2
Long=0.595      #Longuer
h=0.015      #hauteur
h2=h*2      
I=(b*h**3)/12 #Moment d'inertie
I2=(b2*h2**3)/12
EI=E*I       #rigidity
EI2=E*I2
A=h*b        #aire section
A2=b2*h2
m1=Long*b*h*rho #Masse
####


N_cell = 2  #Numero cells
N=int(N_cell/2)
freq = np.logspace(1,3,10000) #frequency
frf=[]
frf_tmd=[]
auto_va_list_2=[]
auto_va_list_3=[]
auto_va_list_2_im_1=[]
auto_va_list_3_im_2=[]
Vr=1
l=Long/N_cell
for f in freq:
    s=2*np.pi*f
    T=[]    
    T2=[]
            
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
    
    #kcoef=1/(1/(k_HSLDS+c2*1j*2*np.pi*f)-1/(4*np.pi**2*f**2*m2))
    
    #F=np.array([[1,0,0,0],[0,1,0,0],[-kcoef,0,1,0],[0,0,0,1]])
    
    kap2=(( (rho*A2) / (EI2) ) * s**2)**(1/4)
    Tii_2=np.cos(l*kap2)+np.cosh(l*kap2)
    T12_2=(np.sin(l*kap2)+np.sinh(l*kap2))/kap2
    T13_2=-(np.sin(l*kap2)-np.sinh(l*kap2))/(EI2*kap2**3)
    T14_2=(np.cos(l*kap2)-np.cosh(l*kap2))/(EI2*kap2**2)
    T21_2=-(np.sin(l*kap2)-np.sinh(l*kap2))*kap2
    T24_2=-(np.sin(l*kap2)+np.sinh(l*kap2))/(EI2*kap2)
    T31_2=(EI2*kap2**3)*(np.sin(l*kap2)+np.sinh(l*kap2))
    T32_2=-(EI2*kap2**2)*(np.cos(l*kap2)-np.cosh(l*kap2))
    T34_2=-T21_2
    T43_2=-T12_2
    T42_2=(EI2*kap2)*(np.sin(l*kap2)-np.sinh(l*kap2))
    
    T2.append(0.5*np.array([[Tii_2,T12_2,T13_2,T14_2],[T21_2,Tii_2,-T14_2,T24_2],[T31_2,T32_2,Tii_2,-T21_2],[-T32_2,T42_2,-T12_2,Tii_2]]))
    
    T_p = T2[-1] @ T[-1]
    
    
    T_p = np.linalg.matrix_power(T_p, N)
    
    
    Vl = Vr/(T_p[2,2] - (T_p[2,3]*T_p[3,2]/T_p[3,3]))        
    Ml = -T_p[3,2]*Vl/(T_p[3,3])
    Yr = (T_p[0,2]*Vl)  +  (T_p[0,3]*Ml)
    

    auto_va, auto_vetores=np.linalg.eig(T_p)

    
    frf.append(abs(Yr))
    f=np.log(frf)
    mu = 1j*np.log(auto_va)/Long
    mu.sort()
    auto_va_list_2.append(mu[2].real)
    auto_va_list_3.append(mu[3].real)
    auto_va_list_2_im_1.append(mu[2].imag)
    auto_va_list_3_im_2.append(mu[3].imag)
    

plt.plot(freq,f,label='Sans absorbeur')
plt.title('Frequency Response Function (FRF) of Euler-Bernoulli Beam')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('FRF Amplitude')
plt.legend()
plt.grid(True)
plt.show()

#print(auto_va_list_2)



plt.plot(freq,auto_va_list_2,".")
plt.plot(freq,auto_va_list_3,".")
plt.title('Auto-valeurs de la matrice de transfert')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Constant de propagation')
plt.grid(True)
plt.show()

plt.plot(freq,auto_va_list_2_im_1,".")
plt.plot(freq,auto_va_list_3_im_2,".")
plt.title('Auto-valeurs de la matrice de transfert')
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Constant de propagation')
plt.grid(True)
plt.show() 

arquivo = 'caminho/do/seu/arquivo.txt'
dados = np.genfromtxt(arquivo, skip_header=1)

# Separar os dados em colunas
freq = dados[:, 0]
log_abs_v = dados[:, 1]

# Plotar o gráfico
plt.plot(freq, log_abs_v, label='Log10(abs(v))')
plt.title('Gráfico Logarítmico')
plt.xlabel('Frequência (Hz)')
plt.ylabel('log10(abs(v))')
plt.grid(True)
plt.legend()
plt.show()