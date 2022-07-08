
% function forvin: definisce i carichi e i vincoli della struttura

function [nVinc,dV,nForze,dF]=forvin;

 %Matrice dei carichi: la i-esima riga di dF contiene il numero del nodo 
 %                     caricato, la direzione del grado di libertà secondo 
 %                     cui agisce la componente della forza e il valore 
 %                     assegnato alla stessa.
 %  dF(i,1)=numero nodo;
 %  dF(i,2)=gdl caricato ("1" per indicare "x" o "2" per indicare "y");
 %  dF(i,3)=valore della forza assegnata;

  dF(1,:)=[2,2,-30000];

  [nForze,nn]=size(dF);  %nForze = numero totale dei carichi applicati



 %Matrice dei vincoli: la i-esima riga di dV contiene il numero del nodo 
 %                     vincolato, la direzione del grado di libertà 
 %                     vincolato e il valore assegnato allo spostamento.
 %  dV(i,1)=numero nodo;
 %  dV(i,2)=gdl vincolato ("1" per indicare "x" o "2" per indicare "y");
 %  dV(i,3)=valore dello spostamento assegnato;

  dV(1,:)=[1,1,0];
  dV(2,:)=[1,2,0];
  dV(3,:)=[3,1,0];
  dV(4,:)=[3,2,0];

  [nVinc,nn]=size(dV);  %nVinc = numero totale dei gdl vincolati

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
