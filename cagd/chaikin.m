function A = chaikin(b, k)
%CHAIKIN za dan nabor tock in parameter k izvede k korakov
%   Chaikinovega postopka

if k == 0
    A = b;
    return
end

[m, n] = size(b);
A = NaN(2*(m-1),n);

iA = 1;
ib = 1;

while ib ~= m
    A(iA,:) = b(ib,:)*3/4 + b(ib+1,:)/4;
    A(iA+1,:) = b(ib,:)/4 + b(ib+1,:)*3/4;
    
    iA = iA+2;
    ib = ib+1;
end

if k ~= 1
    A = chaikin(A, k-1);
end

end

