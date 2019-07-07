function b = polyder3(B, d, T)
% POLYDER3
% Odvod polinoma p dveh spr. stopnje n v smereh d, izvrednoten v tocki T

r = size(d, 1);
n = size(B, 1);

U = NaN(n-1,3);
U(1:r,:) = d;
U(r+1:n-1,:) = repmat(T, [n-r-1, 1]);

b = factorial(n-1) / factorial(n-r-1) * blossom3(B,U);

end

