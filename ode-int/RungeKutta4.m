function [y,x] = RungeKutta4(fun, a, b, y0, h)
% 4. stopenjska Runge-Kutta metoda reda 4 za dif. en.

x = a:h:b;
n = length(x);
y = [y0, zeros(length(y0), n-1)];

for i = 2:n
    k1 = h*fun(x(i-1), y(:,i-1));
    k2 = h*fun(x(i-1) + h/2, y(:,i-1) + k1/2);
    k3 = h*fun(x(i-1) + h/2, y(:,i-1) + k2/2);
    k4 = h*fun(x(i), y(:,i-1) + k3);
    y(:,i) = y(:,i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
end

end


