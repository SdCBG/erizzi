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
%              Laboratorio di giovedi' 2 Maggio 2013                  %
%                                                                     %
%                Trave 2D incastrata ad un estremo                    %
%                                                                     %
%                                             Autore: Luca Gambirasio %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. INTRODUZIONE DEI DATI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A1. INTRODUZIONE DEI DATI GEOMETRICI DELLA STRUTTURA
% INPUT. Numero di nodi
nodi=4;
% INPUT. Numero di travi
travi=3;
% INPUT. Matrice delle coordinate dei nodi
XY=[0,0;
    1000,0;
    2000,0;
    3000,0];
% INPUT. Matrice delle incidenze del sistema
Inc=[1,2;
     2,3;
     3,4];
Inc=[Inc,3*Inc(:,1)-2,3*Inc(:,1)-1,3*Inc(:,1),3*Inc(:,2)-2,3*Inc(:,2)-1,3*Inc(:,2)];
     % Ogni riga è associata ad una trave
     % Prima colonna: Primo nodo della trave
     % Seconda colonna: Secondo nodo della trave
     % Terza colonna: gdl orizzontale del primo nodo della trave
     % Quarta colonna: gdl verticale del primo nodo della trave
     % Quinta colonna: gdl di rotazione del primo nodo della trave
     % Sesta colonna: gdl orizzontale del secondo nodo della trave
     % Settima colonna: gdl verticale del secondo nodo della trave
     % Ottava colonna: gdl di rotazione del secondo nodo della trave

% A2. DEFINIZIONE DEI DATI CINEMATICI DELLA STRUTTURA
% Numero totale di gradi di libertà dei nodi
gdltot=3*nodi;
% INPUT. gdl vincolati (assegnati)
gdlV=[1,2,3]; % E' il nodo a sinistra (incastrato)
% INPUT. Valore degli spostamenti assegnati ai gdl vincolati
UgdlV=[0,0,0];
% gdl liberi
gdlL=1:gdltot;
gdlL(gdlV)=[];

% A3. DEFINIZIONE DELLE CARATTERISTICHE DELLE TRAVI
% INPUT. Vettore delle aree delle travi
A=[100;100;100]; % [mm^2]
% INPUT. Vettore dei moduli di Young delle travi
E=[206000;206000;206000]; % [MPa]
% INPUT. Vettore dei momenti d'inerzia delle travi
I=[10000;10000;10000]; % [mm^4]
% Lunghezze delle travi
L=zeros(travi,1); % Vettore delle lunghezze
for ne=1:travi
    % Nodi dell'ne-esima trave
    nA=Inc(ne,1);   % Primo nodo
    nB=Inc(ne,2);   % Secondo nodo
    % Coordinate dei nodi dell'ne-esima trave
    XA=XY(nA,:);   % Coordinate del primo nodo
    XB=XY(nB,:);   % Coordinate del secondo nodo
    % Lunghezza dell'ne-esima trave
    L(ne,1)=norm(XB-XA);
end

% A4. DEFINIZIONE DEI NODI CARICATI
% Vettore per le possibili forze e coppie applicate ai nodi
forze=zeros(gdltot,1);
% INPUT. Carichi distribuiti applicati alle travi
pq=zeros([travi,2]); % 1° colonna carichi distribuiti assiali, 2° trasversali
pq(:,1)=[0;0;1]; % carichi assiali sulle 3 travi [N/mm]
pq(:,2)=[0.01;0.01;0.01]; % carichi trasversali sulle 3 travi [N/mm]
% Introduzione forze equivalenti a carichi distribuiti nel vettore forze
for ne=1:travi
    fe=zeros([6,1]); % vettore per l'introduzione delle forze equivalenti,
    fe([1,4],1)=[1;1]*pq(ne,1)*L(ne,1)/2; % forze assiali
    fe([2,5],1)=[1;1]*pq(ne,2)*L(ne,1)/2; % forze tangenziali
    fe([3,6],1)=[1;-1]*pq(ne,2)*L(ne,1)^2/12; % coppie
    % Vettore delle incidenze dell'ne-esima trave, dice a quali gdl la
    % trave afferisce
    v=Inc(ne,3:8);
    forze(v,1)=forze(v,1)+fe;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B. RISOLUZIONE DELLA STRUTTURA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% B1. INTRODUZIONE DEL SISTEMA GLOBALE
% Vettore dei termini noti (forze esterne) delle travi
T=forze;
% Definizione del vettore delle RV del sistema
S=zeros(gdltot,1);
% Definizione della matrice di rigidezza del sistema
K=zeros(gdltot,gdltot);
% Calcolo matrici di rigidezza travi e loro introduzione nel sistema
for ne=1:travi
    % Matrice di rigidezza dell'ne-esima trave
    ke=zeros(6,6);
    ke([1,4],[1,4])=[1,-1;-1,1]*E(ne,1)*A(ne,1)/L(ne,1);
    ke([2,5],[2,5])=[12,-12;-12, 12]*E(ne,1)*I(ne,1)/L(ne,1)^3;
    ke([3,6],[3,6])=[4,2;2,4]*E(ne,1)*I(ne,1)/L(ne,1);
    ke([2,5],[3,6])=[6,6;-6,-6]*E(ne,1)*I(ne,1)/L(ne,1)^2;
    ke([3,6],[2,5])=[6,-6;6,-6]*E(ne,1)*I(ne,1)/L(ne,1)^2;
    % Vettore delle incidenze dell'ne-esima trave, dice a quali gdl la
    % trave afferisce
    v=Inc(ne,3:8);
    % Introduzione della matrice di rig. dell'ne-esima trave nel sistema
    K(v,v)=K(v,v)+ke;
end
% Definizione del vettore degli spostamenti e rotazioni del sistema
U=zeros(gdltot,1);
U(gdlV,1)=UgdlV; % necessario se gli spost./rot. assegnati sono diversi da zero

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

% B5. CALCOLO DELLE AZIONI INTERNE
% Vettore delle azioni interne delle travi
Nint=zeros(travi,21);
Tint=zeros(travi,21);
Mint=zeros(travi,21);
for ne=1:travi  
    % Matrice di rigidezza dell'ne-esima trave
    ke=zeros(6,6);
    ke([1,4],[1,4])=[1,-1;-1,1]*E(ne,1)*A(ne,1)/L(ne,1);
    ke([2,5],[2,5])=[12,-12;-12, 12]*E(ne,1)*I(ne,1)/L(ne,1)^3;
    ke([3,6],[3,6])=[4,2;2,4]*E(ne,1)*I(ne,1)/L(ne,1);
    ke([2,5],[3,6])=[6,6;-6,-6]*E(ne,1)*I(ne,1)/L(ne,1)^2;
    ke([3,6],[2,5])=[6,-6;6,-6]*E(ne,1)*I(ne,1)/L(ne,1)^2;   
    % Introduzione forze equivalenti a carichi distribuiti nel vettore forze
    fe=zeros([6,1]);
    fe([1,4],1)=[1;1]*pq(ne,1)*L(ne,1)/2; % forze assiali
    fe([2,5],1)=[1;1]*pq(ne,2)*L(ne,1)/2; % forze tangenziali
    fe([3,6],1)=[1;-1]*pq(ne,2)*L(ne,1)^2/12; % coppie
    % forze di estremità
    v=Inc(ne,3:8); % vettore incidenze, gdl cui la ne-esima trave riferisce
    Ue=U(v,1); % spostamenti dell'ne-esima trave
    R=ke*Ue-fe;
    % Calcolo azione assiale (positiva di trazione) dell'asta ne-esima
    Nint(ne,1)=-R(1,1);
    dL=L(ne,1)/20;
    for np=2:21
        Nint(ne,np)=Nint(ne,np-1)-pq(ne,1)*dL;
    end
    % Calcolo azione di taglio (positiva oraria) dell'asta ne-esima
    Tint(ne,1)=R(2,1);
    for np=2:21
        Tint(ne,np)=Tint(ne,np-1)+pq(ne,2)*dL;
    end
    % Calcolo momento (positivo se tende fibre opposte ad asse y) dell'asta
    % ne-esima
    Mint(ne,1)=-R(3,1);
    for np=2:21
        Mint(ne,np)=Mint(ne,np-1)+Tint(ne,np-1)*dL+pq(ne,2)*(dL^2)/2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C. PRESENTAZIONE DEI RISULTATI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp 'RISOLUZIONE DELLA STRUTTURA:'
disp 'Vettore delle reazioni vincolari del sistema' 
S
disp 'Vettore degli spostamenti e rotazioni del sistema'
U

% C1. RAPPRESENTAZIONE DELLA CONFIGURAZIONE DEFORMATA
figure(1)
title('Configurazione deformata')
hold on
for ne=1:travi
  nA=Inc(ne,1);
  nB=Inc(ne,2);
  XA=XY(nA,1);
  XB=XY(nB,1);
  v=Inc(ne,[3,4,5,6,7,8]);
  Ue=U(v,1);
  dL=L(ne,1)/20;
  ds=[0:.05:1]';
  Udef=[1-ds,ds]*Ue([1,4],1);
  vdef=[1-3*ds.^2+2*ds.^3,L(ne,1)*ds.*(1-2*ds+ds.^2), 3*ds.^2-2*ds.^3, L(ne,1)*ds.*(-ds+ds.^2)]*Ue([2,3,5,6],1);
  Udef=Udef+(4*ds.*(1-ds))*pq(ne,1)*L(ne,1)^2/8/E(ne,1)/A(ne,1);
  vdef=vdef+(16*(ds.^2).*((1-ds).^2))*pq(ne,2)*(L(ne,1)^4)/384/E(ne,1)/I(ne,1);
  plot([XA,XB],[0,0],'k-')
  plot([XA:dL:XB]'+Udef,vdef,'b-')
end

% C2. CREAZIONE DEI DIAGRAMMI DELLE AZIONI INTERNE
% Azione assiale
figure(2)
view([0,0,-1])
title('Azione assiale')
hold on
for ne=1:travi
  dL=L(ne,1)/20;
  nA=Inc(ne,1);
  nB=Inc(ne,2);
  XA=XY(nA,1);
  XB=XY(nB,1);
  fill([XA:dL:XB,XB,XA],[Nint(ne,:),0,0],'r-')
end

% Taglio
figure(3)
view([0,0,-1])
title('Azione di taglio')
hold on
for ne=1:travi
  dL=L(ne,1)/20;
  nA=Inc(ne,1);
  nB=Inc(ne,2);
  XA=XY(nA,1);
  XB=XY(nB,1);
  fill([XA:dL:XB,XB,XA],[Tint(ne,:),0,0],'r-')
end

% Momento flettente
figure(4)
view([0,0,-1])
title('Momento flettente')
hold on
for ne=1:travi
  dL=L(ne,1)/20;
  nA=Inc(ne,1);
  nB=Inc(ne,2);
  XA=XY(nA,1);
  XB=XY(nB,1);
  fill([XA:dL:XB,XB,XA],[Mint(ne,:),0,0],'r-')
end