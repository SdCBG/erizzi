% Corso di Complementi di Scienza delle Costruzioni
% Universita' di Bergamo, Facolta' di Ingegneria, Dalmine
% Docente: prof. Egidio Rizzi
%
% PROGRAMMA PER IL CALCOLO DELLE TENSIONI PRINCIPALI
% scritto da G. Cocchetti e E. Rizzi
% maggio 2006
%
% Inizializzazioni:
% cancella variabili e figure eventualmente in memoria da run precedenti
%
clear all
clf
%
% Dati di input
% Tensore sforzo di Cauchy sigma (simmetrico)
% con valori dati da inserire
%
sig=[40, 10, 20; 
     10, 50, 30; 
     20, 30, 60]
%
% Check per sigma non simmetrico
%
if(norm(sig-sig',1)>0)
    disp('Fornire tensore sforzo simmetrico')
    return
end
%
% Calcolo di sigma^2, sigma^3 e operatori traccia
%
sigq=sig*sig;
sigc=sigq*sig;
trsig=sig(1,1)+sig(2,2)+sig(3,3);
trsigq=sigq(1,1)+sigq(2,2)+sigq(3,3);
trsigc=sigc(1,1)+sigc(2,2)+sigc(3,3);
%
% Calcolo degli invarianti di sigma
%
I1=trsig;
I2=1/2*(trsigq-trsig^2);
I3=1/3*trsigc-1/2*trsig*trsigq+1/6*trsig^3;
%
% Check su I3
%
I3-det(sig);
%
% Deviatore di sforzo e suoi invarianti da calcolo diretto
%
I=eye(3);
p=trsig/3;
sigv=p*I;
s=sig-sigv
sq=s*s;
sc=sq*s;
trsq=sq(1,1)+sq(2,2)+sq(3,3);
trsc=sc(1,1)+sc(2,2)+sc(3,3);
J2=1/2*trsq;
J3=1/3*trsc;
%
% Confronto con relazione tra Ji e Ii
%
J2-(1/3*I1^2+I2);
J3-(1/3*I1*(2/9*I1^2+I2)+I3);
%
% Calcolo della soluzione analitica degli autovalori di s
%
alpha1=1/3*acos(3*sqrt(3)/2*J3/J2^(3/2))
alpha2=alpha1+2/3*pi
alpha3=alpha1-2/3*pi
s1=2*(J2/3)^(1/2)*cos(alpha1)
s2=2*(J2/3)^(1/2)*cos(alpha2)
s3=2*(J2/3)^(1/2)*cos(alpha3)
sig1=s1+p
sig2=s2+p
sig3=s3+p
%
% Ordina gli autovalori in ordine crescente
%
sigord=sort([sig1,sig2,sig3])
%
% Calcolo degli autovettori di sigma e di s
% (si ritengono qui gli autovalori distinti)
%
% Primo autovettore corrispondente ad autovalore sig1
% Si ponga n1=1 e si determinano n2 ed n3 dalle 2 ultime eq.ni
%
n1=1;
n=zeros([3,3]);
for i=1:3
    A=[sig(2,2)-sigord(i), sig(2,3); sig(3,2), sig(3,3)-sigord(i)];
    b=-n1*[sig(1,2); sig(1,3)];
    x=A\b;
    n(:,i)=[1; x]/sqrt(1+x(1)^2+x(2)^2);
end
%
% Ridefinizione del terzo autovettore per ottenere terna principale
% destrorsa
%
if(norm(cross(n(:,1),n(:,2))-n(:,3),1)>0.0001)
    n(:,3)=-n(:,3);
end
if(norm(cross(n(:,1),n(:,2))-n(:,3),1)>0.0001)
    n(:,3)=cross(n(:,1),n(:,2));
end
%
% Mostra l'autovettore
%
n
%
% Verifica con calcolo automatico degli autovalori di sigma e di s
% mediante comando Matlab
%
[vecs,vals]=eig(s)
[vecsig,valsig]=eig(sig)
%
% Fig. 1: Rappresentazione grafica della terna principale
%
figure(1)
clf
hold on
view([10,-5,11])
axis(1*[-1,1,-1,1,-1,1]) 
axis equal
%
% Terna cartesiana di riferimento (in rosso)
%
plot3([0;1],[0;0],[0;0],'r-')
plot3([0;0],[0;1],[0;0],'r-')
plot3([0;0],[0;0],[0;1],'r-')
%
% Terna principale locale (nell'origine, in blu)
%
plot3([0;n(1,1)],[0;n(2,1)],[0;n(3,1)],'b-')
plot3([0;n(1,2)],[0;n(2,2)],[0;n(3,2)],'b-')
plot3([0;n(1,3)],[0;n(2,3)],[0;n(3,3)],'b-')
%
% Vista degli sforzi nello spazio delle tensioni principali
% secondo la decomposizione sig=s+pI
% sig: tensore sforzo     (in verde)
% s:   sforzo deviatorico (in nero)
% pI:  sforzo volumetrico (in blu)
%
% Fig. 2: Vista da prospettiva generica
%
figure(2)
clf
hold on
view([10,-5,5])
axis(130*[-1,1,-1,1,-1,1]) 
axis equal
%
% Assi sig1, sig2, sig3
%
plot3([0;180],[0;0],[0;0],'r-')
plot3([0;0],[0;180],[0;0],'r-')
plot3([0;0],[0;0],[0;180],'r-')
%
% Vettore che rappresenta il tensore sforzo
%
plot3([0;sigord(1)],[0;sigord(2)],[0;sigord(3)],'g-')
%
% Vettore che rappresenta il deviatore su piano deviatorico
% per l'origine e suo traslato su piano deviatorico per pI
%
plot3([0;sigord(1)-p],[0;sigord(2)-p],[0;sigord(3)-p],'k-')
plot3([0+p;sigord(1)-p+p],[0+p;sigord(2)-p+p],[0+p;sigord(3)-p+p],'k:')
%
% Vettore che rappresenta lo sforzo volumetrico pI e suo traslato per s
%
plot3([0;p],[0;p],[0;p],'b-')
plot3([sigord(1)-p;sigord(1)],[sigord(2)-p;sigord(2)],[sigord(3)-p;sigord(3)],'b:')
%
% Piano deviatorico per pI
%
plot3([3*p;0],[0;3*p],[0;0],'c:')
plot3([3*p;0],[0;0],[0;3*p],'c:')
plot3([0;0],[3*p;0],[0;3*p],'c:')
%
% Piano deviatorico per l'origine
%
plot3([2*p;-p],[-p;2*p],[-p;-p],'c-')
plot3([2*p;-p],[-p;-p],[-p;2*p],'c-')
plot3([-p;-p],[2*p;-p],[-p;2*p],'c-')
%
% Fig. 3: Vista di s sul piano deviatorico dall'asse idrostatico
%         (cambia solo il view point)
%
figure(3)
clf
hold on
view([10,10,10])
axis(100*[-1,1,-1,1,-1,1]) 
axis equal
%
% Assi sig1, sig2, sig3
%
plot3([0;180],[0;0],[0;0],'r-')
plot3([0;0],[0;180],[0;0],'r-')
plot3([0;0],[0;0],[0;180],'r-')
%
% Vettore che rappresenta il tensore sforzo
%
plot3([0;sigord(1)],[0;sigord(2)],[0;sigord(3)],'g-')
%
% Vettore che rappresenta il deviatore su piano deviatorico
% per l'origine e suo traslato su piano deviatorico per pI
%
plot3([0;sigord(1)-p],[0;sigord(2)-p],[0;sigord(3)-p],'k-')
plot3([0+p;sigord(1)-p+p],[0+p;sigord(2)-p+p],[0+p;sigord(3)-p+p],'k:')
%
% Vettore che rappresenta lo sforzo volumetrico pI e suo traslato per s
%
plot3([0;p],[0;p],[0;p],'b-')
plot3([sigord(1)-p;sigord(1)],[sigord(2)-p;sigord(2)],[sigord(3)-p;sigord(3)],'b:')
%
% Piano deviatorico per pI
%
plot3([3*p;0],[0;3*p],[0;0],'c:')
plot3([3*p;0],[0;0],[0;3*p],'c:')
plot3([0;0],[3*p;0],[0;3*p],'c:')
%
% Piano deviatorico per l'origine
%
plot3([2*p;-p],[-p;2*p],[-p;-p],'c-')
plot3([2*p;-p],[-p;-p],[-p;2*p],'c-')
plot3([-p;-p],[2*p;-p],[-p;2*p],'c-')
