
% function aziass: calcola l'azione assiale nelle aste

function [dN]=aziass(nAste,nInc,dXY,dPar,du)

%Vettore delle azioni assiali
 dN=zeros([nAste,1]);
  
 for ne=1:nAste
   %Nodi dell'elemento ne-esimo
     n12=nInc(ne,1:2);

   %Coordinate dei nodi dell'elemento ne-esimo
     dXY12=dXY(n12,:);

   %Lunghezza dell'elemento ne-esimo
     dLne=norm(dXY12(2,:)-dXY12(1,:),2);

   %Coseno e seno dell'angolo di inclinazione dell'elemento ne-esimo
     dCSne=(dXY12(2,:)-dXY12(1,:))/dLne;

   %Parametri dell'elemento ne-esimo
     dEAne=dPar(ne,:);

   %Matrice di rigidezza dKne dell'elemento ne-esimo
     [dKne]=matrig(dXY12,dEAne);

   %Vettore dei gradi di libertà dei nodi dell'elemento ne-esimo
     nVne=nInc(ne,3:6);

   %Spostamenti nei nodi dell'elemento ne-esimo
     dUne=du(nVne,1);

   %Forze nei nodi dell'elemento ne-esimo
     dFne=dKne*dUne;

   %Azione assiale dell'elemento ne-esimo
     dN(ne,1)=dCSne*dFne(3:4,1);
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
