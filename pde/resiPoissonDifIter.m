function [U,X,Y,s] = resiPoissonDifIter(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K,U0,T,S,M,w)
% Opis:
%  resiPoissonDifIter z uporabo diferencne metode resi Poissonovo enacbo na
%  pravokotniku pri Dirichletovih robnih pogojih, pri cemer se za racunanje
%  priblizkov uporablja Jacobijeva, Gauss-Seidlova ali SOR metoda.
%
% Definicija:
%  [U,X,Y,s] = resiPoissonDifIter(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K,U0,T,S,M)
%
% Vhodni podatki:
%  a,b,c,d  parametri, ki dolocajo pravokotnik (a,b) x (c,d),
%  F        funkcija, ki doloca Poissonovo enacbo - laplace U = F,
%  Gc,Gd    funkciji, ki dolocata robne pogoje v x smeri:
%           U(x,c) = Gc(x), U(x,d) = Gd(x)
%  Ga,Gb    funkciji, ki dolocata robne pogoje v y smeri:
%           U(a,y) = Ga(y), U(b,y) = Gb(y)
%  J,K      parametra diskretizacije pri diferencni metodi, ki dolocata
%           stevilo notranjih tock mreze v x oziroma y smeri,
%  U0       funkcija, ki doloca zacetne priblizke v notranjosti mreze,
%  T        toleranca maksimalne absolutne razlike med komponentama dveh
%           zaporednih priblizkov v iteraciji,
%  S        maksimalno stevilo korakov iteracije,
%  M        parameter, ki doloca iterativno metodo: 1 = Jacobijeva metoda,
%           2 = Gauss-Seidlova metoda, 3 = SOR metoda (pri optimalnem
%           parametru)
%
% Izhodni podatki:
%  U        tabela velikosti (K+2) x (J+2), ki predstavlja numericno
%           resitev Poissonove enacbe - laplace U = F pri danih robnih
%           pogojih,
%  X,Y      tabeli velikosti (K+2) x (J+2), ki vsebujeta x in y koordinate
%           tock mreze, na kateri se izvede diferencna metoda (vrednost
%           U(k,j) torej predstavlja numericni priblizek za resitev
%           Poissonove enacbe v tocki (X(k,j), Y(k,j)),
%  s        stevilo opravljenih korakov iteracije


% koraka
dx = (b-a)/(J+1);
dy = (d-c)/(K+1);

% oznake
delta_2 = (dx^2 * dy^2) / (2*(dx^2 + dy^2));
% theta_x = delta_2 / dx^2;
% theta_y = delta_2 / dy^2;
theta_x = dy^2 / (2*(dx^2 + dy^2));
theta_y = dx^2 / (2*(dx^2 + dy^2));

% tocke
x = a:dx:b;
y = c:dy:d;
yin = y(2:end-1);

[X,Y] = meshgrid(x,y);
Xin = X(2:end-1,2:end-1);
Yin = Y(2:end-1,2:end-1);

% init
Ur = NaN(K+2,J+2);

Ur(1,:) = Gc(x);
Ur(end,:) = Gd(x);
Ur(2:end-1,1) = Ga(yin');
Ur(2:end-1,end) = Gb(yin');

Ur(2:end-1,2:end-1) = U0(Xin,Yin);

% SOR param.
if M == 3
    if nargin < 16
        rhoJ = 1 - 4 * theta_x * sin(pi/(2*(J+1)))^2 - ...
                   4 * theta_y * sin(pi/(2*(K+1)))^2;

        w = 2 / (1 + sqrt(1 - rhoJ^2));
    end
end

% loop
s = 1;
Ur1 = Ur;

while s < S
    for j = (1:J)+1
        for k = (1:K)+1
            if M == 1
                Ur1(k,j) = theta_x * (Ur(k,j-1) + Ur(k,j+1)) + ...
                           theta_y * (Ur(k-1,j) + Ur(k+1,j)) + ...
                           delta_2 * F(x(j),y(k));
            elseif M >= 2
                Ur1(k,j) = theta_x * (Ur1(k,j-1) + Ur(k,j+1)) + ...
                           theta_y * (Ur1(k-1,j) + Ur(k+1,j)) + ...
                           delta_2 * F(x(j),y(k));
                if M == 3
                    Ur1(k,j) = w*Ur1(k,j) + (1-w)*Ur(k,j);
                end
            end
        end
    end
    
    if max(abs(Ur1 - Ur)) <= T
        break
    else
        Ur = Ur1;
        s = s+1;
    end
end

U = Ur1;

end

