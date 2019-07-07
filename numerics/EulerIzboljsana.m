function [y,x] = EulerIzboljsana(fun, a, b, y0, h)
% Modificirana (dvostopenjska) Eulerjeva eksplicitna metoda za dif. en.

x = a:h:b;
n = length(x);
y = [y0, zeros(length(y0),n-1)];

for i = 2:n
    k1 = h*fun(x(i-1), y(:,i-1));
    k2 = h*fun(x(i-1) + h/2, y(:,i-1) + k1/2);
    y(:,i) = y(:,i-1) + k2;
end

end

