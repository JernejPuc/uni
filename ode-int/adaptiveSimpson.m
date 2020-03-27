function [Q, epsilon] = adaptiveSimpson(f, a, b, delta)
%ADAPTIVESIMSPON
%   [Q, epsilon] = adaptiveSimpson(f, a, b, delta)

%% Integracija za tri tocke oz. n=2
% Dolzina koraka
h = (b-a)/2;

% Vmesna tocka
c = a + (b-a)/2;

% Simpsonovo tretjinsko pravilo
Q1 = h/3*(f(a) + 4*f(c) + f(b));

%% Integracija za pet tock oz. n=4
% Dolzina koraka
h = h/2;

% Novi vmesni tocki
d = a + (c-a)/2;
e = c + (b-c)/2;

% Simpsonovo tretjinsko pravilo
Q2 = h/3*(f(a) + 4*f(d) + 2*f(c) + 4*f(e) + f(b));

%% Adaptivna metoda
% Ocena napake Simpsonovega pravila
epsilon = abs(Q2-Q1)/15;
% epsilon = abs(Q2-Q1)  % Iz navodil za DN

if epsilon > abs(delta)
    % Rekurzivno racunanje podintegralov
    [Qa, epsilona] = adaptiveSimpson(f, a, c, delta/2);
    [Qb, epsilonb] = adaptiveSimpson(f, c, b, delta/2);
    Q = Qa + Qb;
    epsilon = epsilona + epsilonb;
else
    % Boljsi koncni priblizek s korakom Richardsonove ekstrapolacije
    Q = Q2 + (Q2-Q1)/15;
end

end

