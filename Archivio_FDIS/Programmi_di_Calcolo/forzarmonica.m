% Corso di Fondamenti di Dinamica e Instabilita' delle Strutture
% Universita' di Bergamo, Facolta' di Ingegneria, Dalmine
% Docente: prof. Egidio Rizzi
%
% PROGRAMMA PER LA RAPPRESENTAZIONE DELLA RISPOSTA
% DELL'OSCILLATORE SEMPLICE A FORZANTE ARMONICA
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
w1=2*pi/T1;
ke=m*w1^2;
%
% Tempo di analisi (in s) e numero di incrementi temporali
%
ta=10;
ni=1000;
%
% Vettore degli istanti di tempo 
%
Dt=ta/ni;
t=[t0:Dt:ta];
%
% Plot 1 - Rappresentazione del fattore di amplificazione in funzione
%          di beta
%
% Vettore dei valori di beta
%
bvl1=[0:0.001:0.948];
bvg1=[1.050:0.001:3];
Nvl1=abs(1./(1.-bvl1.^2));
Nvg1=abs(1./(1.-bvg1.^2));
%
fig1=figure(1);
set(fig1,'Position',[10 470 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 2]);
plot(bvl1,Nvl1,'b','LineWidth',2)
plot(bvg1,Nvg1,'b','LineWidth',2)
title('Curva di risonanza N(beta) per smorzamento nullo')
xlabel('beta')
ylabel('N(beta)')
%
% Forzante armonica F*sin(w*t)
%
F=1000;
T=ta/5;
w=2*pi/T;
b=w/w1;
%
% Se si vuole assegnare direttamente w rispetto a w1
% (in pratica se si vuole assegnare beta)
%
b=1.2;
w=b*w1;
Ft=F*sin(w*t);
%
% Integrale particolare per beta diverso da 1 
%
ust=F/ke;
U=1/(1-b^2)*ust;
up=U*sin(w*t);
%
% Plot 2 - Rappresentazione dell'integrale particolare
%           per beta<1 in fase con la forzante,
%           per beta>1 in opposizione di fase con la forzante.
%          Nel plot il confronto e' fatto con ust sin(wt) per
%          mostrare l'entita' del fattore di amplificazione
%          e l'andamento della forzante sin(wt), rendendo
%          visibile l'eventuale sfasamento in opposizione di up(t)
%
fig2=figure(2);
set(fig2,'Position',[60 365 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta/10]);
plot(t,up,'b','LineWidth',2)
plot(t,ust*sin(w*t),':r','LineWidth',1.5)
title('Integrale particolare U*sin(wt) riferito a F(t)/k = ust*sin(wt)')
xlabel('t')
ylabel('up(t)')
%
% Integrale particolare in condizioni di risonanza (beta=1)
%
U=-1/2*w1*ust;
upris=U*t.*cos(w1*t);
%
% Plot 3 - Rappresentazione dell'integrale particolare in
%          condizioni di risonanza
%          Nel plot il confronto e' fatto con ust sin(w1t) per
%          mostrare l'ampiezza divergente dell'oscillazione
%          e l'andamento della forzante sin(w1t)
%
fig3=figure(3);
set(fig3,'Position',[110 260 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta/5]);
plot(t,ust*sin(w1*t),':r','LineWidth',1.5)
plot(t,U*t,'--g','LineWidth',1.2)
plot(t,-U*t,'--g','LineWidth',1.2)
plot(t,upris,'b','LineWidth',2)
title('Risonanza: int. part. U*t*cos(w1t) rif. a F(t)/k = ust*sin(w1t)')
xlabel('t')
ylabel('up(t)')
%
% Plot 4 - Illustrazione dei battimenti per beta circa = 1
%
bb=1.1;
uFb=1/(1-bb^2)*ust*(sin((bb*w1)*t)-bb*sin(w1*t));
fig4=figure(4);
set(fig4,'Position',[160 155 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta]);
plot(t,uFb,'b','LineWidth',2)
plot(t,1/(1-bb^2)*ust*2*sin(w1*(bb-1)/2*t),':r','LineWidth',1.5)
plot(t,-1/(1-bb^2)*ust*2*sin(w1*(bb-1)/2*t),':r','LineWidth',1.5)
title('Risposta a forzante per c.i. nulle con battimenti')
xlabel('t')
ylabel('uF(t)')
%
% Integrale totale con risposta alle condizioni iniziali
% Qui sotto riportato un caso con beta < 1
% Variare il parametro beta e osservare la transizione
% battimenti->risonanza per beta -> 1
%
b=1.2;
w=b*w1;
ui=up0/w1*sin(w1*t)+u0*cos(w1*t);
uF=1/(1-b^2)*ust*(sin(w*t)-b*sin(w1*t));
u=ui+uF;
%
% Plot 5 - Rappresentazione della risposta totale a forzante e c.i.
%
fig5=figure(5);
set(fig5,'Position',[210 50 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta/2]);
plot(t,u,'b','LineWidth',2)
plot(t,1/(1-b^2)*ust*sin(w*t),':r','LineWidth',1.5)
title('Risposta totale a F*sin(wt) e c.i. nel caso non smorzato')
xlabel('t')
ylabel('u(t)')
%
%--------------------------------------------------------------------------
% RISPOSTA DELL'OSCILLATORE SMORZATO
%--------------------------------------------------------------------------
%
% Plot 6 - Rappresentazione delle curve di risonanza per vari z
%
% Vettore dei valori di beta
%
bvl1=[0:0.001:0.956];
bvg1=[1.042:0.001:2];
Nvl1=abs(1./(1.-bvl1.^2));
Nvg1=abs(1./(1.-bvg1.^2));
bv=[0:0.001:2];
%z=5%
Nv5=1./((1-bv.^2).^2+(2*(5/100)*bv).^2).^(1/2);
%z=7%
Nv7=1./((1-bv.^2).^2+(2*(7/100)*bv).^2).^(1/2);
%z=10%
Nv10=1./((1-bv.^2).^2+(2*(10/100)*bv).^2).^(1/2);
%z=20%
Nv20=1./((1-bv.^2).^2+(2*(20/100)*bv).^2).^(1/2);
%z=50%
Nv50=1./((1-bv.^2).^2+(2*(50/100)*bv).^2).^(1/2);
%z=1/sqrt(2)=70.71%
Nv70=1./((1-bv.^2).^2+(2*(1/sqrt(2))*bv).^2).^(1/2);
%z=100%
Nv100=1./((1-bv.^2).^2+(2*(100/100)*bv).^2).^(1/2);
%
% Linea traccia dei massimi relativi
%
zvd=[0.05,0.07,0.10,0.20,0.50,0.7071067];
bvmd=sqrt(1-2*zvd.^2);
Nvmd=1./(2*zvd.*(1-zvd.^2).^(1/2));
zv=[0.044:0.01:0.707];
bvm=sqrt(1-2*zv.^2);
Nvm=1./(2*zv.*(1-zv.^2).^(1/2));
%
fig6=figure(6);
set(fig6,'Position',[580 470 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
plot(bvl1,Nvl1,'b','LineWidth',2)
plot(bvg1,Nvg1,'b','LineWidth',2)
plot(bv,Nv5,'r','LineWidth',1)
plot(bv,Nv7,'m','LineWidth',1)
plot(bv,Nv10,'g','LineWidth',1)
plot(bv,Nv20,'c','LineWidth',1)
plot(bv,Nv50,'y','LineWidth',1)
plot(bv,Nv70,'k','LineWidth',1)
plot(bv,Nv100,'b','LineWidth',1)
plot(bvmd,Nvmd,'+','LineWidth',2)
plot(bvm,Nvm,':b','LineWidth',1.1)
title('Curve di risonanza N(beta) per 0 < z < 1')
xlabel('beta')
ylabel('N(beta)')
%
% Plot 7 - Rappresentazione delle curve di fase per vari z
%
bvs=[0:0.001:0.999];
bvd=[1.001:0.001:3];
%z=1%
csivs1=atan(2*(1/100).*bvs./(1-bvs.^2));
csivd1=atan(2*(1/100).*bvd./(1-bvd.^2))+pi;
%z=5%
csivs5=atan(2*(5/100).*bvs./(1-bvs.^2));
csivd5=atan(2*(5/100).*bvd./(1-bvd.^2))+pi;
%z=7%
csivs7=atan(2*(7/100).*bvs./(1-bvs.^2));
csivd7=atan(2*(7/100).*bvd./(1-bvd.^2))+pi;
%z=10%
csivs10=atan(2*(10/100).*bvs./(1-bvs.^2));
csivd10=atan(2*(10/100).*bvd./(1-bvd.^2))+pi;
%z=20%
csivs20=atan(2*(20/100).*bvs./(1-bvs.^2));
csivd20=atan(2*(20/100).*bvd./(1-bvd.^2))+pi;
%z=50%
csivs50=atan(2*(50/100).*bvs./(1-bvs.^2));
csivd50=atan(2*(50/100).*bvd./(1-bvd.^2))+pi;
%z=1/sqrt(2)=70.71%
csivs70=atan(2*(1/sqrt(2)).*bvs./(1-bvs.^2));
csivd70=atan(2*(1/sqrt(2)).*bvd./(1-bvd.^2))+pi;
%z=100%
csivs100=atan(2*(100/100).*bvs./(1-bvs.^2));
csivd100=atan(2*(100/100).*bvd./(1-bvd.^2))+pi;
%
fig7=figure(7);
set(fig7,'Position',[630 365 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
ylim([0 pi]);
set(gca,'YTick',0:pi/4:pi)
set(gca,'YTickLabel',{'0','pi/4','pi/2','3pi/4','pi'})
plot(bvs,0,'b','LineWidth',2)
plot(bvd,pi,'b','LineWidth',2)
plot(bvs,csivs1,'k','LineWidth',1)
plot(bvd,csivd1,'k','LineWidth',1)
plot(bvs,csivs5,'r','LineWidth',1)
plot(bvd,csivd5,'r','LineWidth',1)
plot(bvs,csivs7,'m','LineWidth',1)
plot(bvd,csivd7,'m','LineWidth',1)
plot(bvs,csivs10,'g','LineWidth',1)
plot(bvd,csivd10,'g','LineWidth',1)
plot(bvs,csivs20,'c','LineWidth',1)
plot(bvd,csivd20,'c','LineWidth',1)
plot(bvs,csivs50,'y','LineWidth',1)
plot(bvd,csivd50,'y','LineWidth',1)
plot(bvs,csivs70,'k','LineWidth',1)
plot(bvd,csivd70,'k','LineWidth',1)
plot(bvs,csivs100,'b','LineWidth',1)
plot(bvd,csivd100,'b','LineWidth',1)
title('Curve di fase csi(beta) per 0 < z < 1')
xlabel('beta')
ylabel('csi(beta)')
%
% Integrale particolare nel caso smorzato
%
z=0.05;
b=1.2;
w=b*w1;
N=1/sqrt((1-b^2)^2+(2*z*b)^2);
U=N*ust;
if b<1
 csi=atan(2*z*b/(1-b^2));
elseif b>1
 csi=atan(2*z*b/(1-b^2))+pi;
else b=1 
 csi=pi/2;
end
up=U*sin(w*t-csi);
%
% Plot 8 - Rappresentazione dell'integrale particolare
%          per il sistema smorzato.
%          Nel plot il confronto e' fatto con ust sin(wt) per
%          mostrare l'entita' del fattore di amplificazione
%          e l'andamento della forzante sin(wt), rendendo
%          visibile lo sfasamento rispetto a up(t)
%          (diverso da 0 o Pi nel caso non smorzato)
%
fig8=figure(8);
set(fig8,'Position',[680 260 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta/10]);
plot(t,up,'b','LineWidth',2)
plot(t,ust*sin(w*t),':r','LineWidth',1.5)
title('Integrale part. N*ust*sin(wt-csi) riferito a F(t)/k = ust*sin(wt)')
xlabel('t')
ylabel('up(t)')
%
% Plot 9 - Rappresentazione dell'integrale totale
%          per il sistema smorzato con risposta a forzante e c.i.
%          Si apprezza lo smorzamento della risposta transiente
%          a dare la risposta a regime legata alla forzante armonica.
%          Si noti che anche per beta=1 (risonanza per il sistema non
%          smorzato) le oscillazioni non divergono e si assestano
%          sulla risposta steady-state corrispondente all'integrale
%          particolare
%
up=U*sin(w*t-csi);
wd=w1*sqrt(1-z^2);
D=(1-b^2)^2+(2*z*b)^2;
A=(up0+z*w1*u0+w*ust*(b^2+2*z^2-1)/D)/wd;
B=u0+ust*(2*z*b)/D;
ugoa=exp(-z*w1*t).*(A*sin(wd*t)+B*cos(wd*t));
u=ugoa+up;
%
fig9=figure(9);
set(fig9,'Position',[730 155 560 420]);
axes('XGrid','on','YGrid','on');
box('on');
hold('all');
xlim([0 ta/2]);
plot(t,up,':r','LineWidth',1.5)
plot(t,ugoa,':g','LineWidth',1.5)
plot(t,u,'b','LineWidth',2)
title('Risposta totale a F*sin(wt) e c.i. nel caso smorzato')
xlabel('t')
ylabel('u(t)')
