function [U,X,Y,A,B] = resiPoissonDif(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K)
% Opis:
%  resiPoissonDif z uporabo diferencne metode resi Poissonovo enacbo na
%  pravokotniku pri Dirichletovih robnih pogojih
%
% Definicija:
%  [U,X,Y,A,B] = resiPoissonDif(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K)
%
% Vhodni podatki:
%  a,b,c,d  parametri, ki dolocajo pravokotnik (a,b) x (c,d),
%  F        funkcija, ki doloca Poissonovo enacbo - laplace U = F,
%  Gc,Gd    funkciji, ki dolocata robne pogoje v x smeri:
%           U(x,c) = Gc(x), U(x,d) = Gd(x)
%  Ga,Gb    funkciji, ki dolocata robne pogoje v y smeri:
%           U(a,y) = Ga(y), U(b,y) = Gb(y)
%  J,K      parametra diskretizacije pri diferencni metodi, ki dolocata
%           stevilo notranjih tock mreze v x oziroma y smeri
%
% Izhodni podatki:
%  U        tabela velikosti (K+2) x (J+2), ki predstavlja numericno
%           resitev Poissonove enacbe - laplace U = F pri danih robnih
%           pogojih,
%  X,Y      tabeli velikosti (K+2) x (J+2), ki vsebujeta x in y koordinate
%           tock mreze, na kateri se izvede diferencna metoda (vrednost
%           U(k,j) torej predstavlja numericni priblizek za resitev
%           Poissonove enacbe v tocki (X(k,j), Y(k,j)),
%  A,B      matrika in vektor sistema A*x = B, katerega resitev doloca
%           numericne priblizke v notranjih tockah mreze

% koraka
dx = (b-a)/(J+1);
dy = (d-c)/(K+1);

% oznake
delta_2 = (dx^2 * dy^2) / (2*(dx^2 + dy^2));
% theta_x = delta_2 / dx^2;
% theta_y = delta_2 / dy^2;
theta_x = dy^2 / (2*(dx^2 + dy^2));
theta_y = dx^2 / (2*(dx^2 + dy^2));

% A
if J > 1
    v_x = -theta_x*ones(1,J-1);
    A1 = eye(J) + diag(v_x, -1) + diag(v_x, 1);
    
    A = kron(eye(K), A1);
else
    A = eye(K);
end

if K == 2
    v_y = [zeros(1,J), -theta_y*ones(1,J)];
    A = A + diag(v_y, -J) + diag(fliplr(v_y), J);
elseif K > 2
    v_y = -theta_y*ones(1, J*(K - 1));
    A = A + diag(v_y, -J) + diag(v_y, J);
end

%b
x = a:dx:b;
y = c:dy:d;
xin = x(2:end-1);
yin = y(2:end-1);

xarg = repmat(xin',K,1);
yarg = reshape(repmat(yin,J,1), J*K, 1);
bf = delta_2 * F(xarg, yarg);

if K == 1
    bx = theta_y * (Gc(xin') + Gd(xin'));
elseif K == 2
    bx = theta_y * [Gc(xin'); Gd(xin')];
elseif K > 2
    bx = theta_y * [Gc(xin'); zeros(J*(K-2),1); Gd(xin')];
end

if J == 1
    by = theta_x * (Ga(yin') + Gb(yin'));
elseif J == 2
    by = theta_x * reshape([Ga(yin); Gb(yin)], 1, [])';
elseif J > 2
    by = theta_x * reshape([Ga(yin); zeros(J-2,K); Gb(yin)],1,[])';
end

B = bf + bx + by;

% sis.
sol = A\B;

% out
[X,Y] = meshgrid(x,y);
U = NaN(K+2,J+2);

U(1,:) = Gc(x);
U(end,:) = Gd(x);
U(2:end-1,1) = Ga(yin');
U(2:end-1,end) = Gb(yin');

U(2:end-1,2:end-1) = reshape(sol,J,K)';

end

