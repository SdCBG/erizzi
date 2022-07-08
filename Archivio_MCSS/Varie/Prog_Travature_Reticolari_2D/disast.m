
% function disast: disegna le aste della struttura

function disast(nAste,nInc,dXY,dt,sSt)

 %Disegno delle aste
 for ne=1:nAste
  %Nodi dell'elemento ne-esimo
   n12=nInc(ne,1:2);

  %Coordinate dei nodi dell'elemento ne-esimo
   dx=dXY(n12,1);
   dy=dXY(n12,2);

  %Disegno dell'elemento ne-esimo
   plot(dx,dy,sSt)

   if (dt ~= 0)
    %Lunghezza dell'elemento ne-esimo
     dLne=norm([dx(2)-dx(1),dy(2)-dy(1)],2);

    %Coseno e seno dell'angolo di inclinazione dell'elemento ne-esimo
     dCSne=[dx(2)-dx(1),dy(2)-dy(1)]/dLne;

     dxtext=(dx(1)+dx(2))/2-dt*dCSne(2);
     dytext=(dy(1)+dy(2))/2+dt*dCSne(1);
     text(dxtext,dytext,sprintf('%d',ne),'color',[0,0,1],'linewidth',3)
   end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
