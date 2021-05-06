function [g]=grigio_relativo(C, K, Nmax)
% grigio_relativo determina il vettore g che rende più scuro il colore C
% proporzionalmente al rapporto fra K e Nmax
% 
% Input:    C tripletta di numero fra 0 e 1
%           K quanto si vuole scurire C rispetto a Nmax
%           Nmax massima scurezza: per K= Nmax si ha nero
% Output:   g tripletta di numero fra 0 e 1 tale che C+g è il colore C
%           RGB inscurito proporzionalmente al rapporto fra K ed Nmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=zeros(size(C))-C;
g=K/Nmax*g;
end

