function I32 = I32Op(x, r, Bir, f, s)
%I32OP Summary of this function goes here
%   Detailed explanation goes here

%% Init
I32 = @(x) zeros(size(x));
nx = length(x);
nr = length(r);

% Razsiritev po navodilu: x(-1) = a, x(n+1) = b
X = [x(1), x, x(end)];

% Predvidena simetricnost
I = floor((length(X)-length(x))/2);

%% Operator
for ix = 1:nx
    for ir = 1:nr
        I32 = @(y) I32(y) + (f(x(ix)) + 1/3*s(ix).*(X(ix+I+r(ir))-x(ix))) .* Bir{ix,ir}(y);
    end
end

end

