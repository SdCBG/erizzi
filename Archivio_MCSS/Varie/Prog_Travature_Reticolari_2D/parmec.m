
% function parmec: definisce i parametri meccanici della struttura

function [dPar]=parmec

 %Matrice dei parametri delle bielle: la k-esima riga di dPar contiene 
 %                                    il modulo di Young (E) del materiale 
 %                                    e l'area (A) della sezione 
 %                                    della biella k-esima.
 %  dPar(k,:)=[E, A]

  dE=206000;
  dA=1000;
  dPar=[dE,   dA;
        dE,   dA;
        dE, 3*dA;
        dE,   dA;
        dE,   dA;
        dE, 3*dA;
        dE, 3*dA];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
