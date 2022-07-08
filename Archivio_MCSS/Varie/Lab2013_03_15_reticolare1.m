%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                     Laboratorio del corso di                        %
%                                                                     %
%          Meccanica Computazionale dei Solidi e delle Strutture      %
%                                                                     %
%                                                                     %
%                     prof. Giuseppe COCCHETTI                        %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%             Esercitazione di venerdi' 15 marzo 2013                 %
%                                                                     %
%        Risoluzione di una semplice struttura reticolare             %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     
%                                                                     
%                                                                     
%Definizione dei parametri

E1=206000           %Modulo di Young dell'asta n.1
A1=1000             %Area sezione trasversale dell'asta n.1
L1=2000             %Lunghezza asta n.1
alpha1=0*(pi/180)   %Angolo di inclinazione dell'asta n.1
v1=[3,4,1,2]        %Vettore delle incidenze per l'asta n.1



E2=206000              %Modulo di Young dell'asta n.2
A2=1000                %Area sezione trasversale dell'asta n.2
L2=2000/cos(30/180*pi) %Lunghezza asta n.2
alpha2=30*(pi/180)     %Angolo di inclinazione dell'asta n.2
v2=[5,6,1,2]           %Vettore delle incidenze per l'asta n.2





%Costruzione matrici di rigidezza

%asta n.1
  r1=[cos(alpha1);
      sin(alpha1)]

  k1=E1*A1/L1*[-r1; r1]*[-r1' , r1']


%asta n.2
  r2=[cos(alpha2);
      sin(alpha2)]

  k2=E2*A2/L2*[-r2; r2]*[-r2' , r2']



  
  
%Assemblaggio matrici di rigidezza
K=zeros([6,6])

K(v1,v1)=K(v1,v1)+k1
K(v2,v2)=K(v2,v2)+k2



%Vettore termini noti
T=zeros([6,1])

Fx=30000  %Forza in Newton
Fy=60000  %Forza in Newton

T(1,1)=Fx
T(2,1)=Fy





%Partizione del sistema risolvente
gdl=[1,2]
gdv=[3,4,5,6]

Kuu=K(gdl,gdl)
Kus=K(gdl,gdv)
Ksu=K(gdv,gdl)
Kss=K(gdv,gdv)


Tu=T(gdl,1)
Ts=T(gdv,1)


Us=[0;0;0;0]         %Spostamenti assegnati nei gdv





%Risoluzione del sistema

Uu=Kuu\(Tu-Kus*Us)   %Spostamenti calcolati nei gdl

Ss=Ksu*Uu+Kss*Us-Ts  %Reazioni vincolari nei gdv




%Ri-ordinamento di spostamenti e reazioni vincolari

U=zeros([6,1])  %Spostamenti in tutti i nodi della struttura

U(gdl,1)=Uu     %Spostamenti calcolati nei gdl
U(gdv,1)=Us     %Spostamenti calcolati nei gdv


S=zeros([6,1])  %Reazioni in tutti i nodi della struttura
S(gdv,1)=Ss     %Reazioni calcolate nei gdv





%Calcolo delle azioni assiali nelle aste 

%asta n.1
 N1=E1*A1/L1*[-r1' , r1']*U(v1,1)

%asta n.2
 N2=E2*A2/L2*[-r2' , r2']*U(v2,1)



 
%%%%%  FINE  %%%%
