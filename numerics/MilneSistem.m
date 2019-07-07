function [res,x] = MilneSistem(fun, a, b, y0, h)
% Milne-ova prediktor-korektor metoda reda 4 za dif. en.

x = a:h:b;
n = length(x);

y = RungeKutta4(fun, a, a+3*h, y0, h);
res = [y, zeros(length(y0),n-4)];

for i = 5:n
    fi1 = fun(x(i-1), res(:,i-1));
    fi2 = fun(x(i-2), res(:,i-2));
    fi3 = fun(x(i-3), res(:,i-3));
    yp = res(:,i-4) + 4*h*(2*fi1 - fi2 + 2*fi3)/3;
    res(:,i) = res(:,i-2) + h*(fun(x(i), yp) + 4*fi1 + fi2)/3;
end

end

