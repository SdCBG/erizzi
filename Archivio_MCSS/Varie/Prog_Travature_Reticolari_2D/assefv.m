
% function assefv: assembla carichi e condizioni di vincolo

function [nUs,dUs,nUu,dT]=assefv(nForze,dF,nVinc,dV,nGdlTot)

%Creazione del vettore dei gdl vincolati
% nUs=zeros([nVinc,1]);
% dUs=zeros([nVinc,1]);
% for nv=1:nVinc
%   nn=dV(nv,1);          %Numero del nodo vincolato
%   ni=2*nn-2+dV(nv,2);   %Grado di libertà vincolato
%   nUs(nv,1)=ni;         %Assegnazione del numero del gdl vincolato
%   dUs(nv,1)=dVs(nv,3);  %Assegnazione dell'intensità dello spostamento
% end

  %In alternativa alle precedenti istruzioni: 
   nUs=2*dV(:,1)-2*ones([nVinc,1])+dV(:,2);
   dUs=dV(:,3);


%Riordino dei gdl vincolati e dei corrispondenti spostamenti assegnati
 [nUs,nI]=sort(nUs);
 dUs=dUs(nI,1);


%Definizione del vettore dei gdl non vincolati
 nUu=[1:nGdlTot]';
 nUu(nUs)=[];




%Creazione del vettore dei termini noti (per struttura non vincolata)
 dT=zeros([nGdlTot,1]);

%Assegnazione delle forze 
 for nf=1:nForze
   nn=dF(nf,1);                 %Numero del nodo caricato
   ni=2*(nn-1)+dF(nf,2);        %Grado di libertà caricato
   dT(ni,1)=dT(ni,1)+dF(nf,3);  %Intensità del carico
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
