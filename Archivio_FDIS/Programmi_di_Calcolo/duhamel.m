% Corso di Fondamenti di Dinamica e Instabilita' delle Strutture
% Universita' di Bergamo, Facolta' di Ingegneria, Dalmine
% Docente: prof. Egidio Rizzi
%
% PROGRAMMA PER IL CALCOLO DELL'INTEGRALE DI DUHAMEL
% scritto da G. Cocchetti ed E. Rizzi
% aprile 2006
%
% Inizializzazioni:
% cancella variabili e figure eventualmete in memoria da run precedenti
%
clear all
clf
%
% Dati di input
% Condizioni iniziali a t0 su spostamento u0 e velocita' up0
%
t0=0;
u0=0.5;
up0=1;
%
% Parametri del sistema ad un grado di liberta': massa, fattore
% di smorzamento relativo al critico, rigidezza elastica
%
m=1;
z=15/100;
ke=2000;
w1=sqrt(ke/m);
%
% Tempo di analisi (in s) e numero di incrementi temporali
%
ta=1;
ni=1000;
%
% Warning per ricordare numero pari di intervalli nel
% caso di regola di Simpson
%
if(rem(ni,2)~=0)
 disp('stop: inserire n. di incrementi pari per Simpson)')
 return
end
%
% Vettore degli istanti di tempo
%
Dt=ta/ni;
t=[t0:Dt:ta];
%
% Forzante
%
T=ta/10;
w=2*pi/T;
F=1000;
%
% Forzante armonica
%
Ft=F*cos(w*t);
%
% Forzante generica da inserire
%
% Ft=F*ft(t)
%
% Plot 1 - Rappresentazione della forzante
%
figure(1)
plot(t,Ft)
title('Forzante')
xlabel('t')
ylabel('F(t)')
%
% Vettore degli spostamenti
% Risposta transiente a c.i. non nulle
%
wd=w1*sqrt(1-z^2);
A=(up0+z*w1*u0)/wd;
B=u0;
ut=exp(-z*w1*t).*(A*sin(wd*t)+B*cos(wd*t));
%
% Plot 2 - Rappresentazione della risposta transiente a c.i. non nulle
%
figure(2)
plot(t,ut)
title('Risposta transiente a c.i. non nulle')
xlabel('t')
ylabel('u_t(t)')
%
% Risposta a regime (Integrale di Duhamel)
%
% Regola dei rettangoli
%
ur=zeros(size(ut));
for i=2:(ni+1)
 ti=t(i);
 for k=1:(i-1)
  tk=t(k);
  %
  % funzione h(t) risposta a impulso unitario
  %
  ht=1/m/wd*exp(-z*w1*(ti-tk))*sin(wd*(ti-tk));
  ur(i)=ur(i)+Ft(k)*ht*Dt;
 end
end
%
% Plot 3 - Rappresentazione della risposta a regime (risposta a
% alla forzante secondo Duhamel per c.i. nulle)
%
figure(3)
plot(t,ur)
title('Risposta a regime secondo Duhamel per c.i. nulle')
xlabel('t')
ylabel('u_r(t)')
%
% Plot 4 - Rappresentazione della risposta totale
%
figure(4)
plot(t,ut+ur)
title('Risposta totale')
xlabel('t')
ylabel('u(t)')
%**************************************************************************
% Confronto con risposta analitica
%**************************************************************************
%
% Risposta analitica a forzante armonica F cos(wt)
%
b=w/w1;
D=(1-b^2)^2+(2*z*b)^2;
N=1/sqrt(D);
ust=F/ke;
Z1=ust*(1-b^2)/D;
Z2=ust*(2*z*b)/D;
R=sqrt(Z1^2+Z2^2);
DR=R-ust*N;
csirad=atan(2*z*b/(1-b^2));
csideg=csirad/pi*180;
%
uap=Z1*cos(w*t)+Z2*sin(w*t);
%
% Plot 5 - Rappresentazione della risposta analitica a forzante armonica,
%          integrale particolare in sfasamento sulla forzante
figure(5)
plot(t,uap/Z1,'b')
hold on
plot(t,Ft/F,'g')
title('Risposta analitica a Fcos(wt), integrale particolare con sfasamento sulla forzante')
xlabel('t')
ylabel('u_a_p(t)/u_a_p_0    F(t)/F')
%
% Risposta completa a forzante armonica F cos(wt) con c.i. nulle
% per confronto con l'Integrale di Duhamel calcolato
%
ua0=0;
uap0=0;
Aa=(z*w1*(ua0-Z1)-w*Z2+uap0)/wd;
Ba=ua0-Z1;
uat=exp(-z*w1*t).*(Aa*sin(wd*t)+Ba*cos(wd*t))+uap;
%
% Plot 6 - Rappresentazione della risposta analitica a forzante armonica,
%          integrale totale
figure(6)
plot(t,uat,'b')
hold on
plot(t,ur,'g')
title('Risposta analitica a Fcos(wt) integrale totale, confronto con integrale di Duhamel')
xlabel('t')
ylabel('u_a_t(t)  u_r(t)')