function Bnfab = f2Bnf(f, n, a, b)
%F2BNF Summary of this function goes here
%   Detailed explanation goes here

%% Init. 
if nargin < 4
    a = 0;
end

if nargin < 5
    b = 1;
end

%% Def. Bernstein poly.
Bni = @(n,i,t) nchoosek(n,i) .* t.^i .* (1 - t).^(n-i);

%% Approx.; xn = (b-a)*i/n+a, t = (x-a)/(b-a)
Bnfab = @(x) 0;

for i = 0:n
    Bnfab = @(x) Bnfab(x) + f((b-a)*i/n+a) .* Bni(n,i,(x-a)/(b-a));
end

end

