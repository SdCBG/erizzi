% Corso: Dynamics, Instability and Anelasticity of Structures
% Universita' di Bergamo, Dipartimento di Ingegneria e Scienze Applicate, Dalmine
% Docente: prof. Egidio Rizzi
%
% PROGRAMMA PER L'INTEGRAZIONE DIRETTA DELLE EQUAZIONI DEL MOTO
% TRAMITE IL METODO DI NEWMARK
% Scritto da J. Salvi, E. Rizzi
% Marzo 2015
%
function [u,ud,udd]=newmark_incr(m,c,k,u0,ud0,Ft,dt,t,beta,gamma)
%
% Parametri beta e gamma del metodo - esempi:
% 1) beta=1/4, gamma=1/2: metodo dell'accelerazione media
% (incondizionatamente stabile)
% 2) beta=1/6, gamma=1/2: metodo dell'accelerazione lineare
% (stabile sse dt<0.551*T1)
%
% Condizioni iniziali
%
u=zeros(1,length(t)); ud=u; udd=u;
u(1)=u0; ud(1)=ud0; udd(1)=(Ft(1)-c*ud0-k*u0)/m;
%
% Parametri ausiliari
%
dF=diff(Ft);
keff=k+c*gamma/(beta*dt)+m/(beta*dt^2); 
a=m/(beta*dt)+c*gamma/beta;
b=m/(2*beta)+c*dt*(gamma/(2*beta)-1);
%
% Integrazione delle equazioni del moto
%
for i=1:(length(t)-1);
    dFeff=dF(i)+a*ud(i)+b*udd(i);
    du=dFeff/keff;
    dud=gamma/(beta*dt)*du-gamma/beta*ud(i)+dt*(1-gamma/(2*beta))*udd(i);
    dudd=1/(beta*dt^2)*du-1/(beta*dt)*ud(i)-1/(2*beta)*udd(i);
    u(i+1)=u(i)+du; ud(i+1)=ud(i)+dud; udd(i+1)=udd(i)+dudd;
end