function [a, pf] = ToPoly(p, f, q, x)
%TOPOLY Summary of this function goes here
%   Detailed explanation goes here

%% Init
a = [];

pfi = @(x) zeros(size(x));
pf = {pfi};

%% Loop
for i = 1:length(q)
    qi = q{i};
    
    a = [a, p(f, qi, x)];
    
    pfi = pf{end};
    pfj = @(x) pfi(x) + a(end)*qi(x);
    
    pf{end+1} = pfj;
end

%% Ignore init (zero) placeholder
pf = pf(2:end);

end

