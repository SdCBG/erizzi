
% function disazn: rappresenta il diagramma dell'azione assiale

function disazn(nAste,nInc,dXY,dN,dNmax,dSmax)
 
 for ne=1:nAste
  %Nodi dell'elemento ne-esimo
   n12=nInc(ne,1:2);

  %Coordinate dei nodi dell'elemento ne-esimo
   dx=dXY(n12,1);
   dy=dXY(n12,2);

  %Lunghezza dell'elemento ne-esimo
   dLne=norm([dx(2)-dx(1),dy(2)-dy(1)],2);

  %Coseno e seno dell'angolo di inclinazione dell'elemento ne-esimo
   dCSne=[dx(2)-dx(1),dy(2)-dy(1)]/dLne;

  %Disegno dell'elemento ne-esimo e del relativo diagramma
   dNne=dN(ne);
   if (abs(dN(ne))<dNmax/10^8)
     plot(dx,dy,'g-')
     dNne=0;
   elseif (dN(ne)>0)
     dx=[dx;
         dx(2)-dSmax*abs(dN(ne)/dNmax)*dCSne(2);
         dx(1)-dSmax*abs(dN(ne)/dNmax)*dCSne(2)];
     dy=[dy;
         dy(2)+dSmax*abs(dN(ne)/dNmax)*dCSne(1);
         dy(1)+dSmax*abs(dN(ne)/dNmax)*dCSne(1)];
     fill(dx,dy,'b-')
   else
     dx=[dx;
         dx(2)-dSmax*abs(dN(ne)/dNmax)*dCSne(2);
         dx(1)-dSmax*abs(dN(ne)/dNmax)*dCSne(2)];
     dy=[dy;
         dy(2)+dSmax*abs(dN(ne)/dNmax)*dCSne(1);
         dy(1)+dSmax*abs(dN(ne)/dNmax)*dCSne(1)];
     fill(dx,dy,'r-')
   end
   dxtext=(dx(1)+dx(2))/2+.35*dSmax*dCSne(2);
   dytext=(dy(1)+dy(2))/2-.35*dSmax*dCSne(1);
   text(dxtext,dytext,sprintf('%1.2g',dNne))
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
