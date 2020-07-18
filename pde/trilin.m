function z = trilin(T,t,x,y,o)
% Opis:
%  trilin izracuna vrednosti (odvodov) linearne funkcije, ki je dolocena z
%  vrednostmi v ogliscih trikotnika
%
% Definicija:
%  z = trilin(T,t,x,y,o)
%
% Vhodni podatki:
%  T    tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
%       kartezicni koordinati oglisca trikotnika
%  t    stolpec dolzine 3, v katerem vsak element predstavlja vrednost
%       linearne funkcije v ogliscu trikotnika,
%  x,y  seznama kartezicnih koordinat tock, v katerih racunamo vrednosti
%       funkcije,
%  o    parameter, ki doloca odvod: ce ni podan, metoda vraca vrednosti
%       linearne funkcije, ce je o = 'x' ali o = 'y' pa vrednosti odvodov
%       funkcije po x oziroma po y
%
% Izhodni podatek:
%  z    vrednosti (odvodov) linearne funkcije v tockah, dolocenih s
%       seznamoma x in y

A = [ones(3,1), T];
abc = A \ t;

if nargin < 5
    z = abc(1) + abc(2)*x + abc(3)*y;
elseif o == 'x'
    z = abc(2)*ones(size(x));
elseif o == 'y'
    z = abc(3)*ones(size(y));
end

end
