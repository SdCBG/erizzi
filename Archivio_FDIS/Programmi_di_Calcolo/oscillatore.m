% Corso di Fondamenti di Dinamica e Instabilita' delle Strutture
% Universita' di Bergamo, Facolta' di Ingegneria, Dalmine
% Docente: prof. Egidio Rizzi
%
% PROGRAMMA PER LA RAPPRESENTAZIONE DI OSCILLAZIONI LIBERE
% DELL'OSCILLATORE SEMPLICE
% scritto da E. Rizzi
% aprile 2007
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
u0=1;
up0=20;
%
% Parametri del sistema ad un grado di liberta' non smorzato: 
% massa m, rigidezza elastica Ke
%
m=1;
ke=400; 
w1=sqrt(ke/m);
%
% Volendo assegnare il periodo naturale del sistema:
%
T1=0.5;
Pi=2*asin(1);
w1=2*Pi/T1;
%
% Tempo di analisi (in s) e numero di incrementi temporali
%
ta=1;
ni=1000;
%
% Vettore degli istanti di tempo e degli istanti immediatamente
% successivi a t0 (per rappresentare l'approssimazione u0+u0p*t
% nell'intorno di t0)
%
Dt=ta/ni;
t=[t0:Dt:ta];
ti=[t0:Dt:ta/20];
%
% Oscillazioni libere non smorzate (z=0)
% Vettore degli spostamenti u(t) secondo le varie rappresentazioni
% armoniche alternative
%
A=up0/w1;
B=u0;
R=sqrt(A^2+B^2);
phi=atan(A/B);
psi=atan(B/A);
uz0a=A*sin(w1*t)+B*cos(w1*t);
uz0b=R*cos(w1*t-phi);
uz0c=R*sin(w1*t+psi);
%
% Plot 1 - Rappresentazione delle oscillazione libere non smorzate.
%          Le tre rappresentazioni originano un unico grafico, con
%          sfasamento phi in ritardo rispetto a cos(w1*t) e
%          sfasamento psi in anticipo rispetto a sin(w1*t)
%
fig1=figure(1);
set(fig1,'Position',[225 471 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
plot(t,uz0b,'r','LineWidth',2)
plot(t,uz0c,'g','LineWidth',2)
plot(t,uz0a,'b','LineWidth',2)
plot(t,R*cos(w1*t),':r')
plot(t,R*sin(w1*t),':g')
title('Oscillazioni libere non smorzate (z = 0)')
xlabel('t')
ylabel('u(t)')
%
% Oscillazioni libere smorzate per z<1
% Vettore degli spostamenti u(t) 
%
z=15/100;
wd=w1*sqrt(1-z^2);
A=(up0+z*w1*u0)/wd;
B=u0;
R=sqrt(A^2+B^2);
uzl1=exp(-z*w1*t).*(A*sin(wd*t)+B*cos(wd*t));
%
% Plot 2 - Rappresentazione delle oscillazione libere smorzate per z<1.
%          Si apprezza la differenza col caso non smorzato
%
fig2=figure(2);
set(fig2,'Position',[325 341 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
plot(t,uzl1,'b','LineWidth',2)
plot(t,uz0a,'--r')
plot(t,R*exp(-z*w1*t),':g')
plot(t,-R*exp(-z*w1*t),':g')
title('Oscillazioni libere smorzate per z < 1')
xlabel('t')
ylabel('u(t)')
%
% Risposta smorzata per z=1
% Vettore degli spostamenti u(t) 
%
A=u0;
B=w1*u0+up0;
uz1=exp(-w1*t).*(A+B*t);
%
% Plot 3 - Rappresentazione della risposta per smorz. critico z=1.
%          Si apprezza la differenza col caso non smorzato e col caso
%          con smorzamento inferiore al critico.
%          Viene anche rappresentata la retta u0+up0*t nell'intorno di t0
%
fig3=figure(3);
set(fig3,'Position',[425 211 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
%plot(t,exp(-w1*t),'-.b')
plot(ti,(u0+up0*ti),'-.b')
plot(t,uz1,'b','LineWidth',2)
plot(t,uzl1,':g')
plot(t,uz0a,'--r')
title('Risposta alle c.i. con smorzamento critico z = 1')
xlabel('t')
ylabel('u(t)')
%
% Risposta smorzata per z>1
% Vettore degli spostamenti u(t) 
%
z=2;
lam1=w1*(-z+sqrt(z^2-1));
lam2=w1*(-z-sqrt(z^2-1));
A=-(lam2*u0-up0)/(2*w1*sqrt(z^2-1));
B= (lam1*u0-up0)/(2*w1*sqrt(z^2-1));
uzg1=A*exp(lam1*t)+B*exp(lam2*t);
%
% Plot 4 - Rappresentazione della risposta per smorz. sovracritico z>1.
%          Si apprezza la differenza coi casi precedenti
%
fig4=figure(4);
set(fig4,'Position',[525 81 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
plot(t,uzg1,'b','LineWidth',2)
plot(t,uz1,'-.k')
plot(t,uzl1,':g')
plot(t,uz0a,'--r')
title('Risposta alle c.i. con smorzamento sovracrtico z > 1')
xlabel('t')
ylabel('u(t)')
