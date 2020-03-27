function [res,x] = BDF(fun, a, b, y0, h, varargin)
% s-stopenjska (default 4, max 6) "backward differentiation formula"

if nargin < 6
    s = 4;
else
    s = varargin{1};
end

x = a:h:b;
n = length(x);

y = RungeKutta4(fun, a, a+(s-1)*h, y0, h);
res = [y, zeros(length(y0),n-s)];

opts = optimoptions(@fsolve, 'Display', 'off');

if s == 1
    for i = 2:n
        res(:,i) = fsolve(@(yi) yi - res(:,i-1) - h*fun(x(i), yi), res(:,i-1), opts);
    end
elseif s == 2
    for i = 3:n
        res(:,i) = fsolve(@(yi) yi - 4*res(:,i-1)/3 + res(:,i-2)/3 - 2*h*fun(x(i), yi)/3, res(:,i-1), opts);
    end
elseif s == 3
    for i = 4:n
        res(:,i) = fsolve(@(yi) yi - 18*res(:,i-1)/11 + 9*res(:,i-2)/11 - 2*res(:,i-3)/11 - 6*h*fun(x(i), yi)/11, res(:,i-1), opts);
    end
elseif s == 4
    for i = 5:n
        res(:,i) = fsolve(@(yi) yi - 48*res(:,i-1)/25 + 36*res(:,i-2)/25 - 16*res(:,i-3)/25 + 3*res(:,i-4)/25 - 12*h*fun(x(i), yi)/25, res(:,i-1), opts);
    end
elseif s == 5
    for i = 6:n
        res(:,i) = fsolve(@(yi) yi - 300*res(:,i-1)/137 + 300*res(:,i-2)/137 - 200*res(:,i-3)/137 + 75*res(:,i-4)/137 - 12*res(:,i-5)/137 - 60*h*fun(x(i), yi)/137, res(:,i-1), opts);
    end
else
    for i = 7:n
        res(:,i) = fsolve(@(yi) yi - 360*res(:,i-1)/147 + 450*res(:,i-2)/147 - 400*res(:,i-3)/147 + 225*res(:,i-4)/147 - 72*res(:,i-5)/147 + 10*res(:,i-6)/147 - 60*h*fun(x(i), yi)/147, res(:,i-1), opts);
    end
end

end

