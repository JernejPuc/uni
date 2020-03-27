function [y,x] = EulerEksplicitna(fun, a, b, y0, h)
% Eulerjeva eksplicitna metoda za dif. en.

x = a:h:b;
n = length(x);
y = [y0, zeros(length(y0),n-1)];

for i = 2:n
    y(:,i) = y(:,i-1) + h*fun(x(i-1), y(:,i-1));
end

end

