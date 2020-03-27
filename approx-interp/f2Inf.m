function Infab = f2Inf(f, n, a, b)
%F2INF Summary of this function goes here
%   Detailed explanation goes here

%% Init. 
if nargin < 4
    a = 0;
end

if nargin < 5
    b = 1;
end

%% Def. lin. piece; i = (xn-a)*n/(b-a), xn = (b-a)*i/n+a
h = (b-a)/n;
i = @(x) floor((x-a)/h);

x1 = @(x) h*i(x)+a;
x2 = @(x) h*(i(x)+1)+a;

y1 = @(x) f(x1(x));
y2 = @(x) f(x2(x));

%% Approx.; t = (x-x1)/(x2-x1)
t = @(x) (x-x1(x))./(x2(x)-x1(x));

Infab = @(x) (1-t(x)).*y1(x) + t(x).*y2(x);

end

