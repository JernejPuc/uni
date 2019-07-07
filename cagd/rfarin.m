function q = rfarin(B, w)
% RFARIN - Description
%
% Syntax: q = rfarin(B, w)
%
% Long description

[n, m] = size(B);
q = NaN(n, m);

for i = 1:n-1
    q(i, :) = w(i)/(w(i)+w(i+1))*B(i,:) + w(i+1)/(w(i)+w(i+1))*B(i+1,:);
end
    
end

