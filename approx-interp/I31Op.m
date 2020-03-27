function I31 = I31Op(x, r, Bir, f, df)
%I31OP Summary of this function goes here
%   Detailed explanation goes here

%% Init
I31 = @(x) zeros(size(x));
nx = length(x);
nr = length(r);

% Razsiritev po navodilu: x(-1) = a, x(n+1) = b
X = [x(1), x, x(end)];

% Predvidena simetricnost
I = floor((length(X)-length(x))/2);

%% Operator
for ix = 1:nx
    for ir = 1:nr
        I31 = @(y) I31(y) + (f(x(ix)) + 1/3*df(x(ix)).*(X(ix+I+r(ir))-x(ix))) .* Bir{ix,ir}(y);
    end
end

end

