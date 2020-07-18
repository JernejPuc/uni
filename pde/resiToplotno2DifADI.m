function [U,X,Y] = resiToplotno2DifADI(T,a,b,c,d,alpha,f,ga,gb,gc,gd,N,J,K)
% Opis:
%  resiToplotno2DifADI izracuna resitev toplotne enacbe 
%  du/dt = alpha^2 laplace(u) pri danih zacetnih in robnih pogojih z
%  uporabo diferencne sheme pri predpisani diskretizaciji v kombinaciji z
%  metodo ADI
%
% Definicija:
%  [U,X,Y] = resiToplotno2DifADI(T,a,b,c,d,alpha,f,ga,gb,gc,gd,N,J,K)
%
% Vhodni podatki:
%  T        maksimalni cas,
%  a,b      robni tocki intervala v x smeri, ki doloca pravokotno obmocje
%           (a,b) x (c,d)
%  c,d      robni tocki intervala v y smeri, ki doloca pravokotno obmocje
%           (a,b) x (c,d)
%  alpha    parameter, ki podaja koeficient enacbe,
%  f        funkcija, ki doloca zacetni pogoj,
%  ga,gb    funkciji, ki dolocata robne pogoje vzdolz stranic {a} x (c,d)
%           in {b} x (c,d),
%  gc,gd    funkciji, ki dolocata robne pogoje vzdolz stranic {c} x (a,b)
%           in {d} x (a,b),
%  N        stevilo delilnih tock v casovni smeri (brez 0),
%  J,K      stevilo notranjih delilnih tock v x in y smeri
%
% Izhodni podatki:
%  U        tabela velikosti (K+2)x(J+2)x(N+1), ki vsebuje izracunane
%           priblizke za resitev problema pri casih med 0 in T ter
%           polozajih iz diskretizacije pravokotnika (a,b) x (c,d);
%           natancneje, element tabele U na mestu (k,j,n) predstavlja
%           priblizek pri casu (n-1)*T/N v tocki s kartezicnima
%           koordinatama (a+(j-1)*(b-a)/(J+1), c+(k-1)*(d-c)/(J+2)),
%  X,Y      tabeli kartezicnih koordinat tock, ki dolocajo diskrezizacijo
%           pravokotnika (a,b) x (c,d) in omogocata enostaven prikaz
%           resitve ob izbranem casu (n-1)*T/N z ukazom surf(X,Y,U(:,:,n))

%%
dt = T/N;
dx = (b-a)/(J+1);
dy = (d-c)/(K+1);

t = 0:dt:T;
x_ = a+dx:dx:b-dx;
y_ = c+dy:dy:d-dy;

[X,Y] = meshgrid(linspace(a,b,J+2),linspace(c,d,K+2));
[Xt,Tx] = meshgrid(linspace(a,b,J+2),linspace(0,T,N+1));
[Yt,Ty] = meshgrid(linspace(c,d,K+2),linspace(0,T,N+1));

%%
U = zeros(K+2,J+2,N+1);
U(:,:,1) = f(X,Y);
U(:,1,2:end) = ga(Ty(2:end,:),Yt(2:end,:))';
U(:,end,2:end) = gb(Ty(2:end,:),Yt(2:end,:))';
U(1,:,2:end) = gc(Tx(2:end,:),Xt(2:end,:))';
U(end,:,2:end) = gd(Tx(2:end,:),Xt(2:end,:))';

nux = alpha^2*dt/dx^2;
nuy = alpha^2*dt/dy^2;

%%
Ix = eye(J);
Iy = eye(K);
Nx = diag(2*ones(1,J)) - diag(ones(1,J-1),1) - diag(ones(1,J-1),-1);
Nx = Nx * 0.5 * nux;
Ny = diag(2*ones(1,K)) - diag(ones(1,K-1),1) - diag(ones(1,K-1),-1);
Ny = Ny * 0.5 * nuy;

Gx = @(t) [gc(t,x_); zeros(K-2,J); gd(t,x_)] * 0.5 * nuy;
Gy = @(t) [ga(t,y_)', zeros(K,J-2), gb(t,y_)'] * 0.5 * nux;

%%
for n = 2:N+1
    Un = U(2:end-1,2:end-1,n-1);
    tn = t(n-1);
    Unp12 = ((Ix + Nx)' \ ((Iy - Ny)*Un + Gx(tn) + Gy(tn+dt/2))')';
    Un1 = (Iy + Ny) \ (Unp12*(Ix - Nx) + Gx(tn+dt) + Gy(tn+dt/2));
    U(2:end-1,2:end-1,n) = Un1;
end

end
