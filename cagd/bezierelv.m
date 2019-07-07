function Be = bezierelv(B, k)
% Opis:
%   bezierelv izvede visanje stopnje dane Bezierjeve krivulje
%
% Definicija:
%   Be = bezierelv (B,k)
%
% Vhodna podatka:
%   B   matrika velikosti (n+1) x d, v kateri vsaka vrstica
%       predstavlja d-dimenzionalno kontrolno tocko
%       Bezierjeve krivulje stopnje n,
%   k   stevilo, ki doloca, za koliko zelimo zvisati stopnjo
%       dane Bezierjeve krivulje
%
% Izhodni podatek:
%   Be  matrika velikosti (n+k+1) x d, v kateri vsaka
%       vrstica predstavlja d- dimenzionalno kontrolno tocko
%       Bezierjeve krvulje stopnje n+k, ki ustreza dani
%       Bezierjevi krivulji

% Get dims
[n, d] = size(B);

% Allocate
Be = zeros(n+1, d);

% Interval edges
Be(1, :) = B(1, :);
Be(end, :) = B(end, :);

% Formula for the rest
for i = 2:n
    % Opomba: V formuli je i = 0, 1, ..., n in zato i/(n+1)
    Be(i, :) = (1 - (i-1)/n)*B(i, :) + (i-1)/n*B(i-1, :);
end

% Continue if necessary
if k > 1
    Be = bezierelv(Be, k-1);
end

end

