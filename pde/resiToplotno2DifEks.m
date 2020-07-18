function [U,X,Y] = resiToplotno2DifEks(T,a,b,c,d,alpha,f,ga,gb,gc,gd,N,J,K)
% Opis:
%  resiToplotno2DifEks izracuna resitev toplotne enacbe 
%  du/dt = alpha^2 laplace(u) pri danih zacetnih in robnih pogojih z
%  uporabo eksplicitne diferencne sheme pri predpisani diskretizaciji
%
% Definicija:
%  [U,X,Y] = resiToplotno2DifEks(T,a,b,c,d,alpha,f,ga,gb,gc,gd,N,J,K)
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
for n = 2:N+1
    for j = 2:J+1
        for k = 2:K+1
            U(k,j,n) = U(k,j,n-1) + ...
                       nux*(U(k,j-1,n-1) - 2*U(k,j,n-1) + U(k,j+1,n-1)) + ...
                       nuy*(U(k-1,j,n-1) - 2*U(k,j,n-1) + U(k+1,j,n-1));
        end
    end
end

end
