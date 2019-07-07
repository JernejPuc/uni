function BS = beziersub(B, t, k)
% Opis: beziersub izvede subdivizijo Bezierjeve krivulje
%
% Definicija:
%   BS = beziersub (B,t)
%
% Vhodni podatki :
%   B   matrika kontrolnih tock Bezierjeve krivulje, v
%       kateri vsaka vrstica predstavlja eno kontrolno
%       tocko,
%   t   parameter subdivizije Bezierjeve krivulje
%
% Izhodni podatek :
%   BS  celica, ki vsebuje kontrolne tocke dveh krivulj, ki
%       jih dobimo s subdivizijo prvotne Bezierjeve krivulje

% Get dim
[n, dim] = size(B);

% Alg
D = decasteljau(B, t);

% Allocate partitions
topD = zeros(n, dim);
botD = zeros(n, dim);

% Get partitions from schemes
for i = 1:dim
    topD(:, i) = D{i}(1, :)';
    botD(:, i) = diag(fliplr(D{i}));
end

% Recursion
if k == 1
    BS = {topD botD};
else
    BS = [beziersub(topD, t, k-1) beziersub(botD, t, k-1)];
end

end

