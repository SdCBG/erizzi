%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                       Laboratorio del corso di                      %
%                                                                     %
%          Meccanica Computazionale dei Solidi e delle Strutture      %
%                                                                     %
%                          Giuseppe COCCHETTI                         %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%              Laboratorio di giovedi' 21 marzo 2013                  %
%                                                                     %
%      Codice per la risoluzione di strutture reticolari piane        %
%                                                                     %
%                                             Autore: Luca Gambirasio %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     
%                                                                     
% STRUTTURA DEL CODICE
%   A. INTRODUZIONE DEI DATI
%   B. RISOLUZIONE DELLA STRUTTURA RETICOLARE
%   C. PRESENTAZIONE DEI RISULTATI

% INPUT indica le righe in cui è necessario introdurre dati

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. INTRODUZIONE DEI DATI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A1. INTRODUZIONE DEI DATI GEOMETRICI DELLA STRUTTURA
% Numerare nodi, loro gdl e bielle
% INPUT. Numero di nodi
nodi=5;
% INPUT. Numero di bielle
bielle=7;
% INPUT. Matrice delle coordinate dei nodi
XY=[0,0;
    2000*cos(pi/3),2000*sin(pi/3);
    2000,0;
    2000+2000*cos(pi/3),2000*sin(pi/3);
    2000+2000,0];
% INPUT. Matrice delle incidenze del sistema
Inc=[1,2,1,2,3,4;
     1,3,1,2,5,6;
     2,3,3,4,5,6;
     2,4,3,4,7,8;
     3,4,5,6,7,8;
     3,5,5,6,9,10;
     4,5,7,8,9,10];
     % Ogni riga è associata ad una biella
     % Prima colonna: Primo nodo della biella
     % Seconda colonna: Secondo nodo della biella
     % Terza colonna: gdl orizzontale del primo nodo della biella
     % Quarta colonna: gdl verticale del secondo nodo della biella
     % Quinta colonna: gdl orizzontale del secondo nodo della biella    
     % Sesta colonna: gdl orizzontale del secondo nodo della biella    

% A2. DEFINIZIONE DEI DATI CINEMATICI DELLA STRUTTURA
% Numero totale di gradi di libertà dei nodi
gdltot=2*nodi;
% INPUT. gdl vincolati (assegnati)
gdlV=[1,2,9,10]; % Sono i nodi 1 e 5, in basso alle estremità
% INPUT. Valore degli spostamenti assegnati ai gdl vincolati
UgdlV=[0,0,0,0];
% gdl liberi
gdlL=1:gdltot;
gdlL(gdlV)=[];

% A3. DEFINIZIONE DEI NODI CARICATI
% Vettore per le possibili forze applicate ai nodi
forze=zeros(gdltot,1);
% INPUT. Forze applicate ai gdl dei nodi
forze(6,1)=-50000;   % Forza verticale sul nodo 3, -50000 N

% A4. DEFINIZIONE DELLE CARATTERISTICHE DELLE BIELLE
% INPUT. Vettore delle aree delle bielle
A=[1000;100;100;1000;100;100;1000]; % [mm^2], aste compresse più grosse
% INPUT. Vettore dei moduli di Young delle bielle
E=[206000;206000;206000;206000;206000;206000;206000]; % [MPa]
% Lunghezze e orientamenti delle bielle
L=zeros(bielle,1); % Vettore delle lunghezze
CS=zeros(bielle,2);   % Matrice coseni (1° colonna) e seni (2° colonna)
for ne=1:bielle
    % Nodi dell'ne-esima biella
    nA=Inc(ne,1);   % Primo nodo
    nB=Inc(ne,2);   % Secondo nodo
    % Coordinate dei nodi dell'ne-esima biella
    XA=XY(nA,:);   % Coordinate del primo nodo
    XB=XY(nB,:);   % Coordinate del secondo nodo
    % Lunghezza dell'ne-esima biella
    L(ne,1)=norm(XB-XA);
    % Coseno e seno dell'ne-esima biella
    CS(ne,:)=(XB-XA)/L(ne,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. RISOLUZIONE DELLA STRUTTURA RETICOLARE 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B1. INTRODUZIONE DEL SISTEMA GLOBALE
% Vettore dei termini noti (forze esterne) delle bielle
T=forze;
% Definizione del vettore delle RV del sistema
S=zeros(gdltot,1);
% Definizione della matrice di rigidezza del sistema
K=zeros(gdltot,gdltot);
% Calcolo matrici di rigidezza bielle e loro introduzione nel sistema
for ne=1:bielle
    % Matrice di rigidezza dell'ne-esima biella
    ke=(E(ne,1)*A(ne,1)/L(ne,1))*[CS(ne,:)';-CS(ne,:)']*[CS(ne,:),-CS(ne,:)];
    % Vettore delle incidenze dell'ne-esima biella, dice a quali gdl la
    % biella afferisce
    v=Inc(ne,3:6);
    % Introduzione della matrice di rig. dell'ne-esima biella nel sistema
    K(v,v)=K(v,v)+ke;
end
% Definizione del vettore degli spostamenti del sistema
U=zeros(gdltot,1);
U(gdlV,1)=UgdlV; % necessario se gli spost. assegnati sono diversi da zero

% B2. PARTIZIONE DEL SISTEMA GLOBALE
% Partizione del vettore dei termini noti (forze esterne) del sistema
Tu=T(gdlL,1);
Ts=T(gdlV,1);
% Partizione del vettore delle RV del sistema
Su=S(gdlL,1);
Ss=S(gdlV,1);
% Partizione della matrice di rigidezza del sistema
Kuu=K(gdlL,gdlL);
Kus=K(gdlL,gdlV);
Ksu=K(gdlV,gdlL);
Kss=K(gdlV,gdlV);
% Partizione del vettore degli spostamenti del sistema
Uu=U(gdlL,1);
Us=U(gdlV,1);

% B3. RISOLUZIONE DEL SISTEMA GLOBALE
% Calcolo spostamenti dei gdl liberi
Uu=Kuu\(Tu-Kus*Us);
% Calcolo RV dei gdl vincolati
Ss=Ksu*Uu+Kss*Us-Ts;

% B4. ORDINAMENTO DELLE RV E DEGLI SPOSTAMENTI DEL SISTEMA
% Riordinamento delle RV
S(gdlV,1)=Ss;
% Riordinamento degli spostamenti
U(gdlL,1)=Uu;

% B5. CALCOLO DELLE AZIONI ASSIALI DELLE BIELLE
% Vettore delle azioni assiali nelle bielle
N=zeros(bielle,1);
for ne=1:bielle
    % Vettore delle incidenze dell'ne-esima biella
    v=Inc(ne,3:6);
    % Calcolo dell'azione assiale dell'ne-esima biella
    N(ne,1)=(E(ne,1)*A(ne,1)/L(ne,1))*[-CS(ne,:),CS(ne,:)]*U(v,1); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C. PRESENTAZIONE DEI RISULTATI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'RISOLUZIONE DELLA STRUTTURA RETICOLARE:'
disp 'Vettore delle reazioni vincolari del sistema' 
S
disp 'Vettore degli spostamenti del sistema'
U
disp 'Vettore delle azioni assiali delle bielle del sistema (trazione positiva)'
N




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAPPRESENTAZIONE GRAFICA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Disegno della struttura indeformata
 %Definizione del raggio del cerchio che indica un nodo
  dimXmax=max(XY(:,1))-min(XY(:,1));
  dimYmax=max(XY(:,2))-min(XY(:,2));
  dr=max(dimXmax,dimYmax)/80;
  figure(1)
  clf
  axis([min(XY(:,1))-dimXmax/10,max(XY(:,1))+dimXmax/10,min(XY(:,2))-dimYmax/10,max(XY(:,2))+dimYmax/10]);
  hold on
  axis equal

 %Disegno dei nodi
  for nn=1:nodi
    dxc=XY(nn,1);
    dyc=XY(nn,2);
    nt=[0:pi/20:2*pi];
    dx=dxc+dr*cos(nt);
    dy=dyc+dr*sin(nt);
    plot(dx,dy,'b-')
    text(dxc+2*dr,dyc+2*dr,sprintf('%d',nn),'color',[0,0,1],'linewidth',3)
  end

 %Disegno delle aste
 for ne=1:bielle
   %Nodi dell'elemento ne-esimo
     n1=Inc(ne,1);  %Primo nodo della biella ne-esima
     n2=Inc(ne,2);  %Secondo nodo della biella ne-esima
   %Coordinate dei nodi dell'elemento ne-esimo
     dXn1=XY(n1,1);
     dYn1=XY(n1,2);
     dXn2=XY(n2,1);
     dYn2=XY(n2,2);
     dx=[dXn1,dXn2]';
     dy=[dYn1,dYn2]';
   %Disegno dell'elemento ne-esimo
     plot(dx,dy,'k-')
     text([1,1]*dx/2-2*dr*CS(ne,2),[1,1]*dy/2+2*dr*CS(ne,1),sprintf('%d',ne))
 end


 %Tracciamento del diagramma dell'azione assiale
  drmax=max(dimXmax,dimYmax)/30;
  dNmax=max(abs(N));
  figure(2)
  clf
  axis([min(XY(:,1))-dimXmax/10,max(XY(:,1))+dimXmax/10,min(XY(:,2))-dimYmax/10,max(XY(:,2))+dimYmax/10]);
  hold on
  axis equal
 for ne=1:bielle
   %Nodi dell'elemento ne-esimo
     n1=Inc(ne,1);  %Primo nodo della biella ne-esima
     n2=Inc(ne,2);  %Secondo nodo della biella ne-esima
   %Coordinate dei nodi dell'elemento ne-esimo
     dXn1=XY(n1,1);
     dYn1=XY(n1,2);
     dXn2=XY(n2,1);
     dYn2=XY(n2,2);
     dx=[dXn1,dXn2]';
     dy=[dYn1,dYn2]';
   %Disegno dell'elemento ne-esimo e del relativo diagramma
     if (abs(N(ne))<dNmax/10^8)
       plot(dx,dy,'g-')
     elseif (N(ne)>0)
       dx=[dXn1,dXn2,dXn2-drmax*(abs(N(ne))/dNmax)*CS(ne,2),dXn1-drmax*(abs(N(ne))/dNmax)*CS(ne,2),dXn1]';
       dy=[dYn1,dYn2,dYn2+drmax*(abs(N(ne))/dNmax)*CS(ne,1),dYn1+drmax*(abs(N(ne))/dNmax)*CS(ne,1),dYn1]';
       fill(dx,dy,'b-')
       text([1,1]*dx(1:2,1)/2+drmax/2*CS(ne,2),[1,1]*dy(1:2,1)/2-drmax/2*CS(ne,1),sprintf('%1.2g',N(ne)))
     else
       dx=[dXn1,dXn2,dXn2-drmax*(abs(N(ne))/dNmax)*CS(ne,2),dXn1-drmax*(abs(N(ne))/dNmax)*CS(ne,2),dXn1]';
       dy=[dYn1,dYn2,dYn2+drmax*(abs(N(ne))/dNmax)*CS(ne,1),dYn1+drmax*(abs(N(ne))/dNmax)*CS(ne,1),dYn1]';
       fill(dx,dy,'r-')
       text([1,1]*dx(1:2,1)/2+drmax/2*CS(ne,2),[1,1]*dy(1:2,1)/2-drmax/2*CS(ne,1),sprintf('%1.2g',N(ne)))
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


