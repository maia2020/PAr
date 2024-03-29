close all
clc

num_mode=1;  % N° du mode auquel on s'intéresse 

%Plaque 1 L=59.5cm,b=3cm, h=1.5cm acier

eta_w=0.0003*217*1i;
E=(1+eta_w)*69e9;     %Module de Young
rho=2710;    %masse volumique
b=0.03;      %épaisseur
h=0.015;     %hauteur
I=(b*h^3)/12;%Moment d'inertie
EI=E*I;      %flexural rigidity
A=h*b;       %aire section

ordo=[];      %liste des ordonnées pour la FRF sans absorbeur
ordo_TMD=[];  %liste des ordonnées pour la FRF avec absorbeur

lengthf=10000; %pas de discrétisation pour la fréquence
Vr=1; % amplitude de l'impulstion fixée par défaut à 1 

k_HSLDS=1.1423e+04; %Raideur de l'absorbeur TMD (résultat des différents choix de dimensionnement) 



%__________________________________________________________________________
% Paramètres globaux (à modifier selon la poutre)
L=0.595;%longueur (m)
b=0.03;%largeur (m)
ep=0.015;%épaisseur (m)
rho=2710; %épaisseur (kg/m3)
%m1=L*b*ep*rho;%masse (kg)
%m1=0.44331; %masse modale associée au mode 1
%m1=0.13642; %masse modale assm2ociée au mode 2
%m1=0.047025;%masse modale associée au mode 3

%participation modale
lm1=[0.44331 0.13642 0.047025];
m1=lm1(num_mode);              %plutôt que la masse de la poutre on prend la masse qui participe au mode d'intérêt

%------------fréquence de résonance à isoler------------------
%freq=34.766 % mode 1 
%freq=217.01  % mode 2
%freq=596.4  % mode 3
lfreq=[34.766 217.01 596.4];
freq=lfreq(num_mode);          %fréquence d'intérêt où l'on souhaite appliquer une atténuation 

% f1=34.1598;
% f2=34.9079;
% Q=freq/(f2-f1);
ksi1=0.01;         %amortissement dans la poutre (+/- arbitrairement)
w1=2*pi*freq; 
k1=w1^2*m1;        %raideur de la poutre
c1=2*ksi1*w1*m1;   %coeff d'amortissement

%déformée =>masse modale
%__________________________________________________________________________
% Paramètres dimensionnement
mu=0.05; %ratio de masse m2/m1
ksi2op=sqrt(3*mu/8/(1+mu))+0.1616*ksi1/(1+mu);
qop=1/(1+mu)*(1-1.5906*ksi1*sqrt(mu/(1+mu)));

%==========Caractéristiques du résonateur HSLDS==========
w2=qop*w1;
m2=mu*m1;
c2=2*ksi2op*w2*m2;
k_HSLDS=w2^2*m2;
%========================================================

%========Tracé de la FRF=======
for f=logspace(1,3,lengthf)
        
    
s=2*pi*1i*f;
    %Décomposition de la poutre en 3 segments (une coupure au milieu (0,5 (en pourcentage de la poutre)),
    %une coupure juste juste avant la fin (0,5+0.499+0.999 ; c'est en ce point que sera systématiquement appliquée la force), la dernière au bout
    %de la poutre (0.5+0.499+0.001=1)
    
    %l'intérêt est de pouvoir placer la force où l'on souhaite sur la
    %poutre en changeant les coeffs)

   
    coef=[0.5 0.5];
    l1=L*coef(1);
    l2=L*coef(2);
    l=[l1 l2];
    l_T=[]; %liste des 3 matrices de transfert associée au 3 morceaux de la poutre
    for k=1:1:2
    
    %===========Construction de la matrice de transfert de la poutre seule (pas de TMD)========
   
    kap=((2*pi*f)^2*rho*A/(EI))^(1/4);

    Tii=cos(l(k)*kap)+cosh(l(k)*kap);
    T12=(sin(l(k)*kap)+sinh(l(k)*kap))/kap;
    T13=-(sin(l(k)*kap)-sinh(l(k)*kap))/(EI*kap^3);
    T14=(cos(l(k)*kap)-cosh(l(k)*kap))/(EI*kap^2);
    T21=-(sin(l(k)*kap)-sinh(l(k)*kap))*kap;
    T24=-(sin(l(k)*kap)+sinh(l(k)*kap))/(EI*kap);
    T31=(EI*kap^3)*(sin(l(k)*kap)+sinh(l(k)*kap));
    T32=-(EI*kap^2)*(cos(l(k)*kap)-cosh(l(k)*kap));
    T34=-T21;
    T43=-T12;
    T42=(EI*kap)*(sin(l(k)*kap)-sinh(l(k)*kap));

    Tflex=0.5*[Tii T12 T13 T14
        T21 Tii -T14 T24
        T31 T32 Tii -T21
        -T32 T42 -T12 Tii];
    
    l_T=[l_T Tflex];
    end
    Tl1=l_T(:,1:4); %matrice de transfert du premier morceau (0 à 0.5) 
    Tl2=l_T(:,5:8); %matrice de transfert du 2e morceau (0.5 à 0.999) 
    %Tl3=l_T(:,9:12); %matrice de transfert du dernier morceau (0.999 à 1) 

    kcoef=1/(1/(k_HSLDS+c2*i*2*pi*f)-1/(4*pi^2*f^2*m2)); %raideur perçue par l'extrémité de la poutre avec le TMD au bout de la poutre (voir photo envoyée)
    F=[1 0 0 0
       0 1 0 0
       -kcoef 0 1 0
       0 0 0 1]; %F est l'effet du TMD

    T=Tl2*Tl1; %matrice de transfert globale pour poutre sans TMD
 
    T_TMD=Tl2*F*Tl1;       %matrice de transfert globale pour poutre avec TMD
    
    
    %on extrait maintenant l'amplitude de vibration du bout libre de la poutre
    Vl=Vr/(T(3,3)-T(3,4)*T(4,3)/T(4,4));
    Ml=-T(4,3)/T(4,4)*Vl;
    yr=T(1,3)*Vl+T(1,4)*Ml;
    
    Vl_TMD=Vr/(T_TMD(3,3)-T_TMD(3,4)*T_TMD(4,3)/T_TMD(4,4));
    Ml_TMD=-T_TMD(4,3)/T_TMD(4,4)*Vl_TMD;
    yr_TMD=T_TMD(1,3)*Vl_TMD+T_TMD(1,4)*Ml_TMD;
    
    ordo_TMD=[ordo_TMD abs(yr_TMD)/Vr];
    ordo=[ordo abs(yr)/Vr];  % on pondère par l'amplitude de l'impulsion, ici Vr=1
     
end

figure
%G=ampli/100;
%loglog(fr,G,'linewidth',1.2);
%hold on

x=logspace(1,3,lengthf);
loglog(x,ordo,'linewidth',1.2)
hold on
loglog(x,ordo_TMD,'linewidth',1.2)
hold on
grid on
xlabel('Frequence (Hz)');
ylabel('FRF (dB)');
legend({'Sans résonateur','Avec résonateur'},'FontSize',18)
set(findall(gcf,'type','text'),'FontSize',18,'fontWeight','normal')
