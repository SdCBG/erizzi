
% function matrig: genera la matrice di rigidezza di un'asta

function [dKne]=matrig(dXY12,dEAne);

   %Calcolo della lunghezza dell'elemento ne-esimo
    % dXn1=dXY12(1,1);
    % dYn1=dXY12(1,2);
    % dXn2=dXY12(2,1);
    % dYn2=dXY12(2,2);
    % dLne=sqrt((dXn2-dXn1)^2+(dYn2-dYn1)^2);

   %In alternativa alle precedenti istruzioni: 
     dLne=norm(dXY12(2,:)-dXY12(1,:),2);



   %Coseno e seno dell'angolo di cui è inclinato l'elemento ne-esimo
    % dCne=(dXn2-dXn1)/dLne;
    % dSne=(dYn2-dYn1)/dLne;
    % dCSne=[dCne,dSne];

   %In alternativa alle precedenti istruzioni: 
     dCSne=(dXY12(2,:)-dXY12(1,:))/dLne;



   %Matrice di rigidezza dell'elemento ne-esimo
     dEne=dEAne(1);
     dAne=dEAne(2);

    % dQne=[-dCne,-dSne,dCne,dSne]';
   %In alternativa alla precedente istruzione: 
     dQne=[-dCSne,dCSne]';


     dKne=(dEne*dAne/dLne)*dQne*dQne';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
