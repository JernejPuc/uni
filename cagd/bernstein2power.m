function b = bernstein2power(p)
% Opis:
%   power2bernstein pretvori polinom, predstavljen s
%   koeficienti v potencni bazi, v polinom, predstavljen
%   v Bernsteinovi bazi
%
% Definicija:
%   p = bernstein2power(p)
%
% Vhodni  podatek:
%   b   seznam koeficientov dolzine n+1, ki po vrsti
%       pripadajo razvoju polinoma stopnje n v Bernsteinovi
%       bazi od 0-tega do n-tega Bernsteinovega baznega
%       polinoma
%
% Izhodni  podatek:
%   p   seznam koeficientov dolzine n+1, ki po vrsti
%       pripadajo razvoju polinoma stopnje n v potencni
%       bazi od x^n do 1

% Get poly dim
n = size(p, 2) - 1;

% Reverse order
p = fliplr(p);

% Init out
b = zeros(1, n + 1);

% Apply formula
for i = 0:n
    for j = i:n
       b(i+1) = b(i+1) + (-1)^(i+j) * nchoosek(n, j) * nchoosek(j, i) * p(j+1);
    end
end

end

