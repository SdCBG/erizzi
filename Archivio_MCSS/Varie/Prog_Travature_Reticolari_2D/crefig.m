
% function crefig: crea la finestra per una figura

function crefig(nF,vCoordFig,vCoordAssi,cTit)

  hF=figure(nF);
  clf
  axis(vCoordAssi);
  hold on
  axis equal
  title(cTit);
  set(hF,'Position',vCoordFig);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
