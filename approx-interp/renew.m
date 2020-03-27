function Enew = renew(r, E, y, a, b)
%RENEW Summary of this function goes here
%   Detailed explanation goes here

%% Init.
if nargin < 4
    a = E(1);
end

if nargin < 5
    b = E(end);
end

%% Remez I
if (a <= y) && (y < E(1))
    if sign(r(E(1))) == sign(r(y))
        E(1) = y;
    else
        E(end) = y;
    end
elseif (E(end) <= y) && (y <= b)
    if sign(r(E(end))) == sign(r(y))
        E(end) = y;
    else
        E(1) = y;
    end
else
    xi = find((E-y) <= 0, 1, 'last');
    
    if sign(r(E(xi))) == sign(r(y))
        E(xi) = y;
    else
        E(xi+1) = y;
    end
end

Enew = sort(E);

end

