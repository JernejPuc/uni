function [U,x] = resiPoVetru(c,f,g,h,T,a,b,N,J)
% Opis:
%  resiPoVetru izracuna resitev advekcijske enacbe du/dt + c*du/dx = 0 pri
%  zacetnem pogoju u(0,x) = f(x) ter robnih pogojih u(t,a) = g(t) in
%  u(t,b) = h(t) z uporabo metode po vetru
%
% Definicija:
%  [U,x] = resiPoVetru(c,f,g,h,T,a,b,N,J)
%
% Vhodni podatki:
%  c        funkcija, ki doloca hitrost potovanja vala v odvisnosti od casa
%           in polozaja na x osi,
%  f        funkcija polozaja, ki doloca zacetni pogoj u(0,x) = f(x),
%  g,h      funkciji casa, ki dolocata robna pogoja u(t,a) = g(t) in
%           u(t,b) = h(t),
%  T        maksimalni cas,
%  a,b      parametra, ki omejujeta prostorsko spremenljivko,
%  N        stevilo delilnih tock v casovni smeri (brez 0),
%  J        stevilo notranjih delilnih tock v prostorski smeri
%
% Izhodni podatek:
%  U        tabela, ki vsebuje izracunane priblizke za resitev problema pri
%           casih med 0 in T z razmikom dt ter polozajih med a in b z
%           razmikom dx; natancneje, element na mestu (n,j) predstavlja
%           priblizek pri casu (n-1)*T/N in polozaju a+(j-1)*(b-a)/(J+1),
%  x        seznam delilnih tocke intervala v prostorski smeri, ki omogoca
%           enostaven prikaz resitve ob izbranem casu (n-1)*T/N z ukazom
%           plot(x,U(n,:))

%%
dt = T/N;
dx = (b-a)/(J+1);

x = a:dx:b;
t = 0:dt:T;

%%
U = zeros(N+1,J+2);
U(1,:) = f(x);
U(2:end,1) = g(t(2:end));
U(2:end,end) = h(t(2:end));

%%
for i = 2:N+1
    for j = 2:J+1
        cnj = c(t(i-1),x(j));
        lam = cnj*dt/dx;
        
        if cnj >= 0
            U(i,j) = lam*U(i-1,j-1) + (1-lam)*U(i-1,j);
        else
            U(i,j) = (1+lam)*U(i-1,j) - lam*U(i-1,j+1);
        end
    end
end

end
