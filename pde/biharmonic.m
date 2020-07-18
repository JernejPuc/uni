function bu = biharmonic(u,x,y,h)
% Opis:
%  biharmonic izracuna aproksimacijo za vrednost biharmonicnega operatorja
%  na funkciji, ki je izracunana na podlagi 13-tockovne sheme
%
% Definition:
%  bu = biharmonic(u,x,y,h)
%
% Vhod:
%  u    funkcija dveh spremenljivk,
%  x,y  kartezicni koordinati tocke v domeni,
%  h    razmik aproksimacije v smereh x in y
%
% Izhod:
%  bu   priblizek za vrednost biharmonicnega operatorja na funkciji u,
%       izracunanan v tocki (x,y)

% uxxxx = (u(x-2*h,y) - 4*u(x-h,y) + 6*u(x,y) - 4*u(x+h,y) + u(x+2*h,y))/h^4;
% uyyyy = (u(x,y-2*h) - 4*u(x,y-h) + 6*u(x,y) - 4*u(x,y+h) + u(x+2*h,y))/h^4;
% uxxyy = (u(x-h,y-h) - 2*u(x,y-h) + u(x+h,y-h)...
%          - 2*u(x-h,y) + 4*u(x,y) - 2*u(x+h,y)...
%          + u(x-h,y+h) - 2*u(x,y+h) + u(x+h,y+h))/h^4;
% 
% bu = uxxxx + 2*uxxyy + uyyyy;

bu = (u(x-2*h,y) + u(x+2*h,y) + u(x,y-2*h) + u(x,y+2*h)...
    + 2*(u(x-h,y-h) + u(x-h,y+h) + u(x+h,y-h) + u(x+h,y+h))...
    - 8*(u(x-h,y) + u(x+h,y) + u(x,y-h) + u(x,y+h))...
    + 20*u(x,y))/h^4;

end

