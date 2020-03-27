function [alpha, beta, q] = TriRek(x, k, N, n)
%TRIREK Summary of this function goes here
%   Detailed explanation goes here

%% Init
qm1 = @(x) zeros(size(x));
q0 = @(x) sqrt(1/(k*N)) * ones(size(x));
p = @(g,h,x) k * g(x)' * h(x);

q = {qm1, q0};
alpha = [];
beta = [sqrt(p(@(x) ones(size(x)), @(x) ones(size(x)), x))];

%% Rek. formula
for i = 1:n
    qim2 = q{end-1};
    qim1 = q{end};
    
    qx = @(x) qim1(x) .* x;
    alpha = [alpha, p(qx, qim1, x)];
    
    qi = @(x) (x-alpha(end)).*qim1(x) - beta(end)*qim2(x);
    
    beta = [beta, sqrt(p(qi,qi,x))];
    qi = @(x) qi(x) / beta(end);
    
    q{end+1} = qi;
end

%% Ignore init (zero) placeholder
q = q(2:end);

end

