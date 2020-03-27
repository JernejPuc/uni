function Bir = TriShema(x, iab, ix, r)
%TRISHEMA Summary of this function goes here
%   Detailed explanation goes here

%% Init
% x-i podintervala za trikotno shemo
xa = x(iab);
xb = x(iab+1);

Xab = [xa, xa, xb, xb];

% Razsiritev po navodilu: x(-1) = a, x(n+1) = b
X = [x(1), x, x(end)];

% Predvidena simetricnost
I = floor((length(X)-length(x))/2);

% Pogoji
Birxj = @(i,r,xj) (X(i+I)-X(i+I-r)) / (X(i+I+r)-X(i+I-r)) * Hi(i,x,xj); 
dBirxj = @(i,r,xj) 3 / (X(i+I+r)-X(i+I-r)) * Hi(i,x,xj);

%% Shema
nab = length(Xab);
S = NaN(nab);

% Vrednosti f v interpoliranih tockah
for i = 1:nab
    S(i,1) = Birxj(ix,r,Xab(i));
end

% Preostanek sheme
for i = 1:nab-1
    for j = 1:nab-i
        % Dopolnilo
        if all(Xab(j:j+i) == Xab(j))
            S(j,i+1) = dBirxj(ix,r,Xab(j));
        else
            S(j,i+1) = (S(j+1,i)-S(j,i)) / (Xab(i+j)-Xab(j));
        end
    end
end

% Koef.
A = S(1,:);

%% Poly
baza = {@(x) ones(size(x))};

for i = 1:nab-1
	baza{i+1} = @(x) baza{i}(x) .* (x-Xab(i));
end

Bir = @(x) zeros(size(x));

for i = 1:nab
    Bir = @(x) Bir(x) + A(i)*baza{i}(x);
end

end

