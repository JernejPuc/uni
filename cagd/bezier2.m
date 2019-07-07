function [bx, by, bz] = bezier2(Bx, By, Bz, u, v)
% Opis:
%	bezier2 vrne tocke na Bezierjevi ploskvi iz tenzorskega
%   produkta.
%
% Definicija:
%   [bx, by, bz] = bezier2(Bx, By, Bz, u, v)
%
% Vhodni  podatki:
%   Bx, By, Bz	matrike velikosti n+1 x m+1, ki dolocajo
%               koordinate kontrolnih tock,
%   u, v        vrstici dolzine M in N, ki predstavljata
%               parametre v smereh u in v.
%
% Izhodni  podatki:
%   bx, by, bz  matrike velikosti N x M, ki predstavljajo
%               tocke na Bezierjevi ploskvi:
%               [bx(J,I) by(J,I) bz(J,I)] je tocka pri
%               parametrih u(I)in v(J).

% n should be n+1 but ok
n = size(Bx, 1);
M = length(u);
N = length(v);

q = NaN(n, 3);
bx = NaN(N, M);
by = NaN(N, M);
bz = NaN(N, M);

for i = 1:M
    for j = 1:N
        for k = 1:n
            q(k, :) = bezierder([Bx(k,:)', By(k,:)', Bz(k,:)'], 0, v(j));
        end
        
        p = bezierder(q, 0, u(i));
        bx(i,j) = p(1);
        by(i,j) = p(2);
        bz(i,j) = p(3);
    end
end

end

