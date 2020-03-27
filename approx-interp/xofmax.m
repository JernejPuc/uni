function y = xofmax(r, a, b, N)
%XOFMAX Summary of this function goes here
%   Detailed explanation goes here

%% Init.
if nargin < 4
    N = 1000;
end

%% Mesto max. abs. residuala
x = linspace(a,b,N+1);
rx = r(x);
absrx = abs(rx);
maxrx = max(absrx);

iy = find(absrx == maxrx, 1);
y = x(iy);

end

