
% function geotop: definisce la geometria e la topologia della struttura

function [nInc,nAste,dXY,nNodi]=geotop

%Coordinate dei nodi: la n-esima riga di dXY contiene le coordinate 
% del nodo n-esimo.
% dXY(n,:)=[ascissa del nodo n-esimo,  ordinata del nodo n-esimo]
  dL=2000;
  dXY=[     0,           0;
           dL,           0;
         2*dL,           0;
         dL/2,dL/2*sqrt(3);
       3*dL/2,dL/2*sqrt(3)];

 %Numero totale di nodi
  [nNodi,nn]=size(dXY);  




 %Matrice delle incidenze: la ne-esima riga di nInc contiene i numeri 
 %                         dei nodi posti alle estremità della biella 
 %                         ne-esima e i gradi di libertà associati agli 
 %                         spostamenti dei due nodi.
 %  nInc(ne,:)=[n1, n2, n1u, n1v, n2u, n2v]
 
  nInc=[1, 2;
        3, 2;
        1, 4;
        2, 4;
        2, 5;
        3, 5;
        4, 5];

   nInc=[nInc,nInc(:,1)*2-1,nInc(:,1)*2,nInc(:,2)*2-1,nInc(:,2)*2];


 %Numero totale delle aste
  [nAste,nn]=size(nInc);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
