function u = mke(p,q,r,f,t,g)
% Opis:
%  mke izracuna priblizek za resitev parcialne diferencialne enacbe
%   - d/dx (p(x,y) du/dx) - d/dy (q(x,y) du/dy) + r(x,y) u = f(x,y)
%  z robnim pogojem u = g po metodi koncnih elementov z zveznimi odsekoma
%  linearnimi funkcijami nad triangulacijo
%
% Definicija
%  u = mke(p,q,r,f,t,g)
%
% Vhodni podatki:
%  p,q,r,f  funkcije dveh spremenljivk, ki dolocajo parcialno diferencialno
%           enacbo,
%  t        triangulacija obmocja, predstavljena z razredom triangulation,
%  g        funkcija dveh spremenljivk, ki doloca vrednost resitve na robu
%           obmocja
%
% Izhodni podatek:
%  u        3D triangulacija, predstavljena z razredom triangulation, ki
%           doloca zvezno odsekoma linearno funkcijo, ki je priblizek za
%           resitev robnega problema po metodi najmanjsih kvadratov

%% Init
% stevilo robnih stranic in tock
bE = freeBoundary(t);
bV = unique(bE(:));

n = size(t.Points,1);    % stevilo tock triangulacije
m = n-length(bV);        % stevilo notranjih tock triangulacije

% vrednosti v tockah triangulacije
U = zeros(n,1);

% indeksiranje notranjih tock
%  ce je v(j) > 0, je j indeks notranje tocke,
%  ce je v(j) = 0, je j indeks robne tocke
v = zeros(n,1);
i = 1;

for j = 1:n
    if ismember(j,bV)
        % vrednosti na robu
        U(j) = g(t.Points(j,1), t.Points(j,2));
    else
        v(j) = i;
        i = i+1;
    end
end

%% MKE
TRI = t.ConnectivityList;
A = zeros(m,m);
b = zeros(m,1);

for tr = 1:size(TRI,1)
    % indeksi tock v trikotniku tr
    tri = TRI(tr,:);
    
    % koordinate tock v trikotniku tr
    trc = [t.Points(tri(1),:); t.Points(tri(2),:); t.Points(tri(3),:)];
    
    for i = 1:3
        % ce notranja tocka
        if ~ismember(tri(i),bV)
            
            % delta v notranji tocki
            ti = zeros(3,1);
            ti(i) = 1;
            
            % indeks glede na A in b
            vi = v(tri(i));
            
            % del notr. tocki pripadajoce bazne funkcije na trikotniku tr
            phii = @(x,y) trilin(trc, ti, x, y);
            dphiidx = trilin(trc, ti, 0, 0, 'x');
            dphiidy = trilin(trc, ti, 0, 0, 'y');
            
            % integrandi za i,i
            int1 = @(x,y) p(x,y)*dphiidx^2 + q(x,y)*dphiidy^2 + r(x,y).*phii(x,y).^2;
            int2 = @(x,y) f(x,y).*phii(x,y);
            
            % prispevek integracije po trikotniku tr
            A(vi,vi) = A(vi,vi) + triintegral(int1, trc);
            b(vi) = b(vi) + triintegral(int2, trc);
            
            for j = i+1:3
                % ce sosedna tocka prav tako notranja
                if ~ismember(tri(j),bV)
                    
                    % delta v sosedni tocki
                    tj = zeros(3,1);
                    tj(j) = 1;
                    
                    % indeks glede na A in b
                    vj = v(tri(j));
                    
                    % del tocki pripadajoce bazne funkcije na trikotniku tr
                    phij = @(x,y) trilin(trc,tj,x,y);
                    dphijdx = trilin(trc,tj,0,0,'x');
                    dphijdy = trilin(trc,tj,0,0,'y');
                    
                    % integrand za i,j
                    int1 = @(x,y) p(x,y)*dphiidx*dphijdx + q(x,y)*dphiidy*dphijdy + r(x,y).*phii(x,y).*phij(x,y);
                    
                    % prispevek integracije po trikotniku tr
                    int1eval = triintegral(int1, trc);
                    
                    A(vi,vj) = A(vi,vj) + int1eval;
                    A(vj,vi) = A(vj,vi) + int1eval;
                end
            end
             
            for j = [1:i-1, i+1:3]
                % ce sosedna tocka na robu
                if ismember(tri(j),bV)
                    
                    % delta v robni tocki
                    tj = zeros(3,1);
                    tj(j) = g(trc(j,1), trc(j,2));
                    
                    % del tocki pripadajoce bazne funkcije na trikotniku tr
                    phi0 = @(x,y) trilin(trc,tj,x,y);
                    dphi0dx = trilin(trc,tj,0,0,'x');
                    dphi0dy = trilin(trc,tj,0,0,'y');
                    
                    % prispevek integracije po trikotniku tr
                    int1 = @(x,y) p(x,y)*dphiidx*dphi0dx + q(x,y)*dphiidy*dphi0dy + r(x,y).*phii(x,y).*phi0(x,y);
                    
                    b(vi) = b(vi) - triintegral(int1, trc);
                end
            end
        end
    end
end

% koeficienti baznih funkcij
alpha = A \ b;

%% Izhodna triangulacija
X = t.Points(:,1);
Y = t.Points(:,2);
U(v~=0) = alpha;

u = triangulation(t.ConnectivityList, X, Y, U);

end
