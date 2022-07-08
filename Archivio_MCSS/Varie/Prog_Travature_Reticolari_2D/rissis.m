
% function rissis: risolve il sistema assemblato di equilibrio

function [du,dR]=rissis(dK,dT,nUu,nUs,dUs,nGdlTot)

%Definizione delle sotto-matrici del sistema
 %Partizione della matrice di rigidezza dK
  dKuu=dK(nUu,nUu);
  dKus=dK(nUu,nUs);
  dKsu=dK(nUs,nUu);
  dKss=dK(nUs,nUs);

 %Partizione del vettore termini noti dT
  dTu=dT(nUu,1);
  dTs=dT(nUs,1);




%Risoluzione del sistema
 %Calcolo spostamenti 
  dUu=dKuu\(dTu-dKus*dUs);

 %Calcolo reazioni vincolari 
  dRs=dKsu*dUu+dKss*dUs-dTs;




%Ricostruzione del vettore degli spostamenti secondo la numerazione dei gdl 
% della struttura
  du=zeros([nGdlTot,1]);
  du(nUu,1)=dUu;
  du(nUs,1)=dUs;

%Ricostruzione del vettore delle reazioni vincolari secondo la numerazione 
% dei gdl della struttura
  dR=zeros([nGdlTot,1]);
  dR(nUs,1)=dRs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
