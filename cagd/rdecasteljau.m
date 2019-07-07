function b = rdecasteljau(B, w, t)
% Opis:
%   rdecasteljau vrne tocko na racionalni Bezierjevi
%   krivulji, ki je izracunana z de Casteljaujevim
%   postopkom, prirejenim za racionalni primer
%
% Definicija:
%   b = rdecasteljau(B,w,t)
%
% Vhodni podatki:
%   B   matrika, katere vrstica predstavlja koordinate
%       kontrolne tocke racionalne Bezierjeve krivulje,
%   w   seznam utezi racionalne  Bezierjeve krivulje,
%   t   stevilo, ki doloca vrednost delilnega parametra v
%       de Casteljaujevem postopku
%
% Izhodni podatek:
%   b   vrstica, ki predstavlja tocko na racionalni
%       Bezierjevi krivulji pri parametru t

% Get dims
[m, d] = size(B);

% Allocate
W = NaN(m);
D = NaN(m, m, d);

% Init
W(1, :) = w;
D(1, :, :) = B;

% Alg
for r = 2:m
    for i = 1:m-r+1
        W(r, i) = (1-t)*W(r-1, i) + t*W(r-1, i+1);
        D(r, i, :) = (1-t)*(W(r-1, i)/W(r, i)) * D(r-1, i,:) + t*(W(r-1, i+1)/W(r, i))*D(r-1 ,i+1, :); 
    end
end

% Endpoints
b = reshape(D(m, 1, :), [1, d]);

end

