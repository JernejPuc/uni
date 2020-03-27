function P = PBaza(x, r)
%PBAZA Summary of this function goes here
%   Detailed explanation goes here

%% Init
nx = length(x);
nr = length(r);
P = cell(nx*nr, nx-1);

%% Loop
for ix = 1:nx
    for ir = 1:nr
        for iab = 1:nx-1
            ixr = (ix-1)*nr + (ir-1) + 1;
            P{ixr, iab} = TriShema(x, iab, ix, r(ir));
        end
    end
end

end

