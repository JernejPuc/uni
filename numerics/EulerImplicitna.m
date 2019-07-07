function [y,x] = EulerImplicitna(fun, a, b, y0, h)
% Eulerjeva implicitna metoda za dif. en.

x = a:h:b;
n = length(x);
y = [y0, zeros(length(y0),n-1)];
opts = optimoptions(@fsolve, 'Display', 'off');

for i = 2:n
    y(:,i) = fsolve(@(yi) yi - y(:,i-1) - h*fun(x(i), yi), y(:,i-1), opts);
end

end

