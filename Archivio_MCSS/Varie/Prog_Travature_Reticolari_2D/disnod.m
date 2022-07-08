
% function disnod: disegna i nodi della struttura

function disnod(nNodi,dXY,dr)

 %Disegno dei nodi
  for nn=1:nNodi
    dxc=dXY(nn,1);
    dyc=dXY(nn,2);
    nt=[0:pi/20:2*pi];
    dx=dxc+dr*cos(nt);
    dy=dyc+dr*sin(nt);
    plot(dx,dy,'r-')
    text(dxc+2*dr,dyc+2*dr,sprintf('%d',nn),'color',[1,0,0],'linewidth',3)
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
