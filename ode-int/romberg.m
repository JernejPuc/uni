function I = romberg(f, a, b, n)

I = zeros(n);

for i=1:n
    I(i,1) = intSimpson13(f, a, b, 2^(i-1));    
end

for j = 2:n
   for i = j:n
       I(i,j) = (4^(j-1) * I(i,j-1) - I(i-1,j-1) ) / (4^(j-1) - 1);
   end
end


end

