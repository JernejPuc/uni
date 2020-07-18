function [y,x] = resiRobniProblem(a,b,p,q,r,alpha,beta,N)
% Opis:
%  resiRobniProblem uporabi diferencno metodo za izracun numericne resitve
%  navadne diferencialne enacbe -(p(x) y')' + q(x) y = r(x) na intervalu
%  (a,b) pri robnih pogojih y(a) = alpha in y(b) = beta
%
% Definicija:
%  [y,x] = resiRobniProblem(a,b,p,q,r,alpha,beta,N)
%
% Vhodni podatki:
%  a, b         robova intervala, na katerem resujemo enacbo,
%  p, q, r      funkcije, ki nastopajo v enacbi,
%  alpha, beta  robni vrednosti resitve,
%  N            parameter, ki doloca ekvidistantno delitev intervala [a,b]
%               na N+1 podintervalov
%
% Izhodna podatka:
%  y            seznam dolzine N+2, ki predstavlja priblizek za resitev
%               robnega problema v N+1 ekvidistantno razporejenih tockah na
%               intervalu [a,b],
%  x            seznam dolzine N+2, ki predstavlja ekvidistantno
%               razporejene tocke na intervlau [a,b]

% centralne tocke
h = (b-a)/(N+1);
x = a:h:b;
xin = x(2:end-1);

% A
diag_lower = -p(xin(2:end) - h/2);
diag_central = p(xin-h/2) + q(xin)*h^2 + p(xin+h/2);
diag_upper = -p(xin(1:end-1) + h/2);

A = diag(diag_lower,-1) + diag(diag_central,0) + diag(diag_upper,1);

% b
b = h^2 * r(xin);
b(1) = b(1) + p(xin(1) - h/2)*alpha;
b(end) = b(end) + p(xin(end) + h/2)*beta;

% sis.
y = (A\b')';
y = [alpha, y, beta];

end

