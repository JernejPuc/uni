function [TRI,U] = trimeshgrid(d)
% Opis:
%  trimeshhgrid sestavi triangulacijo trikotnika z oglisci (0,0), (1,0),
%  (0,1) in izracuna tocke triangulacije v obliki baricentricnih koordinat
%
% Definicija:
%  [TRI,U] = trimeshgrid(d)
%
% Vhodni podatek:
%  d    stevilo intervalov, na katere je razdeljena stranica trikotnika
%
% Izhodna podatka:
%  TRI  matrika velikosti d^2 x 3, v kateri vsaka vrstica predstavlja en
%       trikotnik enakomerne triangulacije, ki deli osnovni trikotnik,
%  U    matrika velikosti (d+2)*(d+1)/2 x 3, v kateri vsaka vrstica
%       predstavlja baricentricne koordinate ene tocke triangulacije

TRI = zeros(d^2,3);
U = zeros((d+1)*(d+2)/2,3);

k = 1;
l = 1;
for i = 0:d
    for j = 0:d-i
        U(k,:) = [d-i-j j i]/d;
        if i < d && j < d-i
            shift = d+1-i;
            TRI(l,:) = [k k+1 k+shift];
            l = l+1;
            if j < d-i-1
                TRI(l,:) = [k+1 k+shift k+1+shift];
                l = l+1;
            end
        end
        k = k+1;
    end
end

end

