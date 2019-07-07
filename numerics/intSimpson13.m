function I = intSimpson13(f, a, b, n)

h = (b - a)/n;
x = linspace(a, b, n+1);
y = f(x);

w = ones(1, n+1);
w(2:2:end-2) = 4;
w(3:2:end-1) = 2;

I = h/3 * dot(w,y);

end

