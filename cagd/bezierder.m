function db = bezierder (B, r, t)
% Opis :
%   bezierder vrne tocke na krivulji , ki predstavlja odvod
%   dane Bezierjeve krivulje
%
% Definicija :
%   db = bezierder (B,r,t)
%
% Vhodni podatki :
%   B   matrika kontrolnih tock Bezierjeve krivulje , v
%       kateri vsaka vrstica predstavlja eno kontrolno
%       tocko ,
%   r   stopnja odvoda , ki ga racunamo ,
%   t   seznam parameterov , pri katerih racunamo odvod
%
%   Izhodni podatek :
%   db  matrika , v kateri vsaka vrstica predstavlja tocko
%       r- tega odvoda pri istoležnem parametru iz seznama t

% Deg. of undiff. poly
n = size(B, 1) - 1;

% No. of coords
d = size(B, 2);

% No. of out pts
k = size(t, 2);

% Init out
db = zeros(k, d);

% Apply formula for every point and dim
for i = 1:k
    D = decasteljau(B, t(i));
    
    for j = 1:d
        % 0th order case (original curve)
        if r == 0
            db(i, j) = D{j}(1, end);
        else
            db(i, j) = factorial(n)/factorial(n-r) * diff(D{j}(1:r+1, n+1-r), r);
        end 
    end
end

end

