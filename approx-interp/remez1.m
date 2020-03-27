function g = remez1(f, a, b, basis, k, n, N, tol, E)
%REMEZ1 Summary of this function goes here
%   Detailed explanation goes here

%% Init.
if nargin < 6
    n = k-1;
end

if nargin < 7
    N = 1000;
end

if nargin < 8
    tol = 1e-6;
end

if nargin < 9
    E = linspace(a,b,n+2);
end

%% Iter.
fprintf('\n');

for i = 1:k
    % p koef.
    [g,m] = minimaxsol(f, basis, E);
    
    p = @(x) (basis(x') * g)';
    
    % Mesto max. abs. err.
    r = @(x) f(x) - p(x);
    y = xofmax(r, a, b, N);
    
    % Report max. err. in diff. glede na minimax
    ae = abs(r(y));
    am = abs(m);
    
    fprintf('i=%d\ne: %f\nd: %f\n', i, ae, ae-am);
    
    % Term crit.
    if abs(r(y))-m < tol
        break
    end
    
    % Zamenjava
    E = renew(r, E, y, a, b);

end

%% Alter.
figure(2);
hold on;
x = linspace(a,b,N+1);
plot(x,r(x));
plot(E,r(E),'o');
plot(x,zeros(size(x)),'--');
title('2d: residual');
hold off;

end

