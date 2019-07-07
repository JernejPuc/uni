function [y,x] = Heunova(fun, a, b, y0, h)
% Heunova metoda za dif. en.

x = a:h:b;
n = length(x);
y = [y0, zeros(length(y0),n-1)];

for i = 2:n
    k1 = h*fun(x(i-1), y(:,i-1));
    k2 = h*fun(x(i), y(:,i-1) + k1);
    y(:,i) = y(:,i-1) + (k1 + k2)/2;
end

end

