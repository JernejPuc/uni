function [d, D] = dividedDifferences(x, y)

n = size(x,2);
d = zeros(n);
d(:,1) = y;

for i = 1:n-1
    for j = 1:n-i
       d(j,i+1) = (d(j+1, i) - d(j, i)) / (x(i+j) - x(j));
    end
end

D = d(1,end);

end

