function [Bd, we] = rbezierelv(B, w)
% Opis:
%   rbezierelv  izvede  visanje  stopnje  dane  racionalne
%   Bezierjeve  krivulje
%
% Definicija:
%   [Be, we] = rbezierelv(B, w)
%
% Vhodna podatka:
%   B   matrika  velikosti (n+1) x d, v kateri  vsaka  vrstica
%       predstavlja d-dimenzionalno  kontrolno  tocko
%       racionalne  Bezierjeve  krivulje  stopnje n,
%   w   seznam  utezi  racionalne  Bezierjeve  krivulje
%
% Izhodni  podatek:
%   Be  matrika  velikosti n+2 x d, v kateri  vsaka  vrstica
%       predstavlja d-dimenzionalno  kontrolno  tocko
%       racionalne  Bezierjeve  krivulje  stopnje n+1, ki je
%       prirejena  dani  racionalni  Bezierjevi  krivulji ,
%   we  seznam  dolzine n+2, v katerem  vsak  element
%       predstavlja  utezi  racionalne  Bezierjeve  krvulje
%       stopnje n+1, ki je  prirejena  dani  racionalni
%       Bezierjevi krivulji

% Get dims
[n, d] = size(B);

% Allocate
Bp = NaN(n, d+1);

% Poly with weighted ctrl pts of (wi*Bi, wi) form, space R^(d+1)
Bp(:, 1:d) = diag(w)*B;
Bp(:, d+1) = w(:);

% Elevate
Be = bezierelv(Bp, 1);

% Get new weights for rational curve from (wi*Bi, wi) form
we = Be(:, d+1);

% Divide ctrl pts back into pts for rational curve (Bi/wi), space R^d
Bd = diag(1./we)*Be(:, 1:d);

end

