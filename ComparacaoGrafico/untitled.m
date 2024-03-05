close all
clc

eta_w=0.0003*217;
E=(1+eta_w)*2.1*e9;     %Module de Young
rho=7850;    %masse volumique
b=0.01;      %épaisseur
h=0.001;     %hauteur
I=(b*h^3)/12;%Moment d'inertie
EI=E*I;      %flexural rigidity
A=h*b;       %aire section
L=0.15 %length

ordo=[];      %liste des ordonnées pour la FRF sans absorbeur
ordo_TMD=[];  %liste des ordonnées pour la FRF avec absorbeur

lengthf=10000; %pas de discrétisation pour la fréquence
Vr=1; % amplitude de l'impulstion fixée par défaut à 1 
kt=4.2

%========Tracé de la FRF=======
for f=logspace(1,3,lengthf)
        
    
s=2*pi*f;
    coef=[0.5 0.499 0.001];
    l1=L*coef(1);
    l2=L*coef(2);
    l3=L*coef(3);
    l=[l1 l2 l3];
    l_T=[]; %liste des 3 matrices de transfert associée au 3 morceaux de la poutre
    for k=1:1:3
    
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
    Tl3=l_T(:,9:12); %matrice de transfert du dernier morceau (0.999 à 1) 

    kcoef=1/(1/(k_HSLDS+c2*i*2*pi*f)-1/(4*pi^2*f^2*m2)); %raideur perçue par l'extrémité de la poutre avec le TMD au bout de la poutre (voir photo envoyée)
    F=[1 0 0 0
       0 1 0 0
       -kcoef 0 1 0
       0 0 0 1]; %F est l'effet du TMD

    T=Tl3*Tl2*Tl1 %matrice de transfert globale pour poutre sans TMD
    break
    T_TMD=Tl3*F*Tl2*Tl1;       %matrice de transfert globale pour poutre avec TMD
    
    
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