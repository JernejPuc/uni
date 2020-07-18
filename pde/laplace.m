function lu = laplace(u,x,y,h)
% Opis:
%  laplace izracuna aproksimacijo za vrednost Laplaceovega operatorja na
%  funkciji, ki je izracunana na podlagi 5-tockovne sheme
%
% Definicija:
%  lu = laplace(u,x,y,h)
%
% Vhod:
%  u    funkcija dveh spremenljivk,
%  x,y  kartezicni koordinati tocke v domeni,
%  h    razmik aproksimacije v smereh x in y
%
% Izhod:
%  lu   priblizek za vrednost Laplaceovega operatorja na funkciji u,
%       izracunanan v tocki (x,y)

% uxx = (u(x-h,y) - 2*u(x,y) + u(x+h,y))/h^2;
% uyy = (u(x,y-h) - 2*u(x,y) + u(x,y+h))/h^2;
% 
% lu = uxx + uyy;

lu = (u(x-h,y) - 4*u(x,y) + u(x+h,y) + u(x,y-h) + u(x,y+h))/h^2;

end

