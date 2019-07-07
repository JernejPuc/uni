function b = decasteljau3(B,U)
% Opis:
%	decasteljau3  izracuna  vrednost  polinoma  dveh  spremenljivk pri
%   parametrih U  = (u,v w)
%
% Definicija:
%   b = decasteljau3(B,U)
%
% Vhodna  podatka:
%   B   matrika  velikosti n+1 x n+1, ki  predstavlja
%       koeficiente  polinoma  dveh  spremenljivk  stopnje n v
%       Bezierjevi  obliki (element  matrike  na  mestu (i,j),
%       j  <= n+2-i, doloca  koeficient  polinoma z indeksom
%       (n+2-i-j, j-1, i-1)),
%   U vektor velikosti 1 x 3, v kateri  vrstice
%       predstavljavo  baricentricne  koordinate  tock  glede
%       na  domenski  trikotnik , za  katere  izvajamo  razcvet
%       polinoma
%
% Izhodni  podatek:
%   b   vrednost polinoma , dolocenega z matriko B,
%       v tockah , dolocenih z vektorjem u

b = blossom3(B, repmat(U, [size(B,1)-1, 1]));

end

