function b = blossom3(B, U)
% Opis:
%   blossom3  izracuna  razcvet  polinoma  dveh  spremenljivk
%
% Definicija:
%   b = blossom3(B, U)
%
% Vhodna  podatka:
%   B   matrika  velikosti n+1 x n+1, ki  predstavlja
%       koeficiente  polinoma dveh  spremenljivk  stopnje n v
%       Bezierjevi  obliki (element  matrike  na  mestu (i,j),
%       j <= n+2-i, doloca  koeficient  polinoma z indeksom
%       (n+2-i-j, j-1, i-1)),
%   U   matrika  velikosti n x 3, v kateri vrstice
%       predstavljavo baricentricne koordinate tock glede
%       na domenski trikotnik, za  katere  izvajamo  razcvet
%       polinoma
%
% Izhodni  podatek:
%   b   vrednost razcveta  polinoma, dolocenega z matriko B,
%       v tockah, dolocenih z matriko u

n = size(B, 1);

% r0
b = B;

% Oziroma r = 2:n, ker je r0 v bistvu prvi primer
for r = 1:n-1
    u = U(r,1);
    v = U(r,2);
    w = U(r,3);
   
    % Zgornje trikotna matrika - r0 ze imamo, torej do n-1
    for i = 1:n-1
        % Meja zopet zaradi zgornje trikotne matrike
        % Po prevedbi (i,j) -> (i',j',k') se vedno velja i+j+k=n-r
        for j = 1:n-i
            b(i,j) = u*b(i,j) + v*b(i,j+1) + w*b(i+1,j);
        end
    end
end

% V resnici primer bn(0,0,0)
b = b(1,1);

end

