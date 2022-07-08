%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              retic2DSt.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                Analisi statica di strutture reticolari                  %
%                                                                         %
%                            Giuseppe COCCHETTI                           %
%                                                                         %
%                          versione del 20/03/2007                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Il presente software è prodotto ed è da utilizzarsi per finalità      %
%   di tipo esclusivamente didattico.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Il codice è costituito dal presente nucleo principale, il quale sfrutta % 
% le chiamate delle seguenti functions:                                   %
%  - assefv: assembla carichi e condizioni di vincolo;                    %
%  - aziass: calcola l'azione assiale nelle aste;                         %
%  - crefig: crea la finestra per una figura;                             %
%  - disast: disegna le aste della struttura;                             %
%  - disazn: rappresenta il diagramma dell'azione assiale;                %
%  - disnod: disegna i nodi della struttura;                              %
%  - forvin: definisce i carichi e i vincoli della struttura;             %
%  - geotop: definisce la geometria e la topologia della struttura;       %
%  - matrig: genera la matrice di rigidezza di un'asta;                   %
%  - parmec: definisce i parametri meccanici della struttura;             %
%  - rissis: risolve il sistema assemblato di equilibrio;                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Analisi statica di strutture reticolari bi-dimensionali')
clear all
format short e


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: acquisizione dei parametri strutturali %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Acquisizione dei parametri di tipo geometrico e topologico
[nInc,nAste,dXY,nNodi]=geotop;

%Acquisizione dei parametri di tipo meccanico
[dPar]=parmec;

%Acquisizione dei carichi e dei vincoli
[nVinc,dV,nForze,dF]=forvin;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elaborazione dei dati del sistema strutturale %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Valutazione del numero totale dei gradi di libertà
% della struttura non vincolata
 nGdlTot=2*nNodi;

%Creazione della matrice di rigidezza (della struttura non vincolata)
 dK=zeros([nGdlTot,nGdlTot]); 

%Ciclo sugli elementi (per assemblaggio matrici di rigidezza dei singoli elementi)
 for ne=1:nAste
   %Numeri dei nodi dell'elemento ne-esimo
    % n1=nInc(ne,1);  %Primo nodo della biella ne-esima
    % n2=nInc(ne,2);  %Secondo nodo della biella ne-esima
    % n12=[n1,n2];
   %In alternativa alle precedenti istruzioni: 
     n12=nInc(ne,1:2);

   %Coordinate dei nodi dell'elemento ne-esimo
    % dXn1=dXY(n1,1);
    % dYn1=dXY(n1,2);
    % dXn2=dXY(n2,1);
    % dYn2=dXY(n2,2);
    % dXY12=[dXn1,dYn1; 
    %        dXn2,dYn2];
   %In alternativa alle precedenti istruzioni: 
     dXY12=dXY(n12,:);

   %Parametri dell'elemento ne-esimo
    % dEne=dPar(ne,1);  %Modulo di Young del materiale dell'elemento ne-esimo
    % dAne=dPar(ne,2);  %Area della sezione trasversale dell'elemento ne-esimo
    % dEAne=[dEne,dAne];
   %In alternativa alle precedenti istruzioni: 
     dEAne=dPar(ne,:);

   %Costruzione della matrice di rigidezza dKne dell'elemento ne-esimo
     [dKne]=matrig(dXY12,dEAne);

   %Assemblaggio della matrice di rigidezza dell'elemento ne-esimo
     nVne=nInc(ne,3:6);                  %Vettore dei gradi di libertà dei nodi dell'elemento ne-esimo
     dK(nVne,nVne)=dK(nVne,nVne)+dKne;   %Operazione di "assemblaggio" della matrice di rigidezza locale dKne
 end

%Creazione del vettore dei vincoli e del vettore dei termini noti (assegnazione dei carichi) 
 [nUs,dUs,nUu,dT]=assefv(nForze,dF,nVinc,dV,nGdlTot);


 
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: generazione dei risultati %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Risoluzione del sistema   
 [du,dR]=rissis(dK,dT,nUu,nUs,dUs,nGdlTot);

%Calcolo delle azioni assiali delle aste 
 [dN]=aziass(nAste,nInc,dXY,dPar,du);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: rappresentazione grafica dei risultati %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Definizione delle dimensioni caratteristiche della struttura
  dXmin=min(dXY(:,1));
  dXmax=max(dXY(:,1));
  dYmin=min(dXY(:,2));
  dYmax=max(dXY(:,2));
  dimXmax=dXmax-dXmin;
  dimYmax=dYmax-dYmin;


%Disegno della struttura indeformata
 %Creazione finestra
  %nF=numero figura
  %vCoordFig =  ascissa e ordinata (in pixel) dell'angolo in basso a sinistra e 
  %             ascissa e ordinata (in pixel) dell'angolo in alto a destra della finestra
  %vCoordAssi = ascissa e ordinata dell'angolo in basso a sinistra e 
  %             ascissa e ordinata dell'angolo in alto a destra della figura
  vCoordFig=[94,517,560,420]; 
  vCoordAssi=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10];
  cTit='Schema strutturale';
  crefig(1,vCoordFig,vCoordAssi,cTit);

 %Definizione del raggio del cerchio che rappresenta un nodo
  dr=max(dimXmax,dimYmax)/80;

 %Disegno dei nodi
  disnod(nNodi,dXY,dr);

 %Definizione della distanza tra etichetta e asta
  dt=max(dimXmax,dimYmax)/40;

 %Disegno delle aste
  disast(nAste,nInc,dXY,dt,'b-');





%Disegno della struttura deformata sovrapposta a quella indeformata
 %Creazione finestra
  vCoordFig=[665,13,560,420]; 
  vCoordAssi=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10]; 
  cTit='Deformata elastica';
  crefig(2,vCoordFig,vCoordAssi,cTit);

 %Calcolo fattore di amplificazione per rappresentazione grafica
  dSmax=max(dimXmax,dimYmax)/30;
  dUmax=max(abs(du));
  dAmplif=10^ceil(log10(dSmax/dUmax));
  text((dXmin+dXmax)/2,dYmax+dimYmax/3,sprintf('fatt. amplif. spost.: %0.5g',dAmplif))

 %Calcolo nuove coordinate dei nodi (con spostamenti amplificati)
  dXYd=zeros(size(dXY));
  for nn=1:nNodi
    dXYd(nn,1)=dXY(nn,1)+du(2*nn-1,1)*dAmplif;
    dXYd(nn,2)=dXY(nn,2)+du(2*nn,1)*dAmplif;
  end

 %Disegno configurazione indeformata
  disast(nAste,nInc,dXY,0,'k-');
  disast(nAste,nInc,dXYd,0,'g--');





%Tracciamento del diagramma dell'azione assiale
 %Creazione finestra
  vCoordFig=[665,517,560,420]; 
  vCoordAssi=[dXmin-dimXmax/10,dXmax+dimXmax/10,dYmin-dimYmax/10,dYmax+dimYmax/10]; 
  cTit='Azione assiale';
  crefig(3,vCoordFig,vCoordAssi,cTit);

 %Calcolo del fattore di amplificazione per la rappresentazione grafica
  dNmax=max(abs(dN));
  dSmax=max(dimXmax,dimYmax)/20;

 %Tracciamento del diagramma
  disazn(nAste,nInc,dXY,dN,dNmax,dSmax);
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

