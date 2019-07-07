function D = decasteljau(b, t)
% Opis:
%   decasteljau  vrne  shemo  de  Casteljaujevega  postopka  za dan
%   seznam  koordinat b pri  danem  parametru t
%
% Definicija:
%   D = decasteljau(b,t)
%
% Vhodna  podatka:
%   b   seznam koordinat kontrolnih tock Bezierjeve krivulje
%       stopnje n,
%   t   parameter, pri  katerem racunamo koordinato
%       Bezierjeve krivulje
%
% Izhodni  podatek:
%   D   tabela  velikosti n+1 x n+1, ki  predstavlja  de
%       Casteljaujevo shemo za koordinate b pri parametru t
%       (element na mestu (1,n+1) je  koordinata Bezierjeve
%       krivulje pri parametru t, elementi na  mestih (i,j)
%       za i > n-j+2 so NaN)

% Get dims
[m, d] = size(b);

% Curve degree
n = m-1;

% Init output
D = cell(1, d);

% Per dim
for i = 1:d
    % Allocate
    D{i} = NaN(m);
    
    % Init
    for j = 0:n
        D{i}(j+1, 1) = b(j+1, i);
    end
    
    % Alg
    for r = 1:n
        for j = 0:n-r
            D{i}(j+1, r+1) = (1-t) * D{i}(j+1, r) + t * D{i}(j+2, r);
        end
    end 
end

end

