function [p,m] = minimaxsol(f, basis, E)
%MINIMAXSOL Summary of this function goes here
%   Detailed explanation goes here

%% Matricni sistem lin. enacb; p(xi) + m*(-1)^i = f(xi)
A = [basis(E'), (-1).^(0:length(E)-1)'];
b = f(E');

res = A\b;
p = res(1:end-1);
m = res(end);

end

