function Hixj = Hi(i,x,xj)
%HI Summary of this function goes here
%   Detailed explanation goes here

%% Init
nx = length(x);
Hixj = 0;

%% Conds.

if i > 1 && x(i-1) <= xj && xj < x(i)
    Hixj = (xj-x(i-1)) / (x(i)-x(i-1));

elseif i < nx && x(i) <= xj && xj < x(i+1)
    Hixj = (x(i+1)-xj) / (x(i+1)-x(i));
    
elseif i == 1 && x(i) <= xj && xj < x(i+1)
    Hixj = (x(i+1)-xj) / (x(i+1)-x(i));
     
elseif i == nx && x(i-1) <= xj && xj <= x(i)
    Hixj = (xj-x(i-1)) / (x(i)-x(i-1));
    
end

end

