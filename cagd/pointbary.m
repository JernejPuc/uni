function [u, v, w] = pointbary(T1, T2, T3, T)
% Pretvori tocko T(x,y) v baricentricne koordinate,
% ki so dolocene z oglisci trikotnika T1, T2, in T3

A = [1 1 1;
    T1(1), T2(1), T3(1);
    T1(2), T2(2), T3(2)];

b = [1; T(1); T(2)];
sol = A \ b;

u = sol(1);
v = sol(2);
w = sol(3);

end

