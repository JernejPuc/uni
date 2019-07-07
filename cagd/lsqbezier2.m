function [Bx, By, Bz] = lsqbezier2(m, n, P, u, v)
% Opis:
%   lsqbezier2 vrne kontrolne to?ke Bezierjeve ploskve, ki
%   se po metodi najmanjaih kvadratov najbolje prilega danim
%   podatkom
%
% Definicija:
%   [Bx,By,Bz] = lsqbezier2(m,n,P,u,v)
%
% Vhodni podatki:
%   m, n        parametra, ki dolocata stopnji Bezierjeve
%               ploskve iz tenzorskega produkta,
%   P           matrika, ki v vrsticah vsebuje tocke v prostoru,       
%               za katere zelimo, da jih Bezierjeva ploskev
%               cim bolje aproksimira,
%   u, v        seznama parametrov, ki dolocata, kateri tocki v
%               domeni pripada posamezna tocka iz P
%
% Izhodni podatki:
%   Bx, By, Bz  matrike velikosti n+1 x m+1, ki predstavljajo
%               kontrolne tocke Bezierjeve ploskve iz
%               tenzorskega produkta, ki se po metodi najmanjaih
%               kvadratov najbolje prilega podatkom

% Stevilo tock
l = size(P, 1);

% M bo imel za vsako tocko (m+1)*(n+1) mesanih produktov baznih polinomov,
% na podlagi katerih bo generiran sistem Ax = b oz. MB = P
% Iscemo take tocke B, za katere bo vsota produktov cim blizje P
% Produkti izhajajo iz izraza:
% [Px, Py, Pz] = Sum(i->m) Sum(j->n) [ B_ij * B(u_i)*B(v_j) ]
M = NaN(l, (m+1)*(n+1));

for k = 1:l
    for j = 0:m
        % Bernsteinov polinom stopnje m, param. z u
        Bu = nchoosek(m,j) * u(k)^j * (1-u(k))^(m-j);
        
        for i = 0:n
            % Bernsteinov polinom stopnje n, param. z v
            Bv = nchoosek(n,i) * v(k)^i * (1-v(k))^(n-i);
            
            M(k, i*(m+1) + j+1) = Bu * Bv;
        end
    end
end

% Resevanje predolocenega sistema MB = P po komponentah
x = M \ P(:,1);
Bx = reshape(x, m+1, n+1)';

y = M \ P(:,2);
By = reshape(y, m+1, n+1)';

z = M \ P(:,3);
Bz = reshape(z, m+1, n+1)';

end

