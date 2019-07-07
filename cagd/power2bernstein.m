function p = power2bernstein(b)
% Opis:
%   power2bernstein  pretvori  polinom , predstavljen s
%   koeficienti v potencni  bazi , v polinom , predstavljen
%   v Bernsteinovi  bazi
%
% Definicija:
%   b = power2bernstein(p)
%
% Vhodni  podatek:
%   p   seznam  koeficientov  dolzine n+1, ki po  vrsti
%       pripadajo  razvoju  polinoma  stopnje n v potencni
%       bazi od x^n do 1
%
% Izhodni  podatek:
%   b   seznam  koeficientov  dolzine n+1, ki po  vrsti
%       pripadajo  razvoju  polinoma  stopnje n v Bernsteinovi
%       bazi od 0-tega do n-tega  Bernsteinovega  baznega
%       polinoma

% Get poly dim
n = size(b, 2) - 1;

% Init output
p = zeros(1, n + 1);

% Apply formula
for i = 0:n
    for j = i:n
       p(i+1) = p(i+1) + nchoosek(j, i) / nchoosek(n, i) * b(j+1);
    end
end

% Reverse order
p = fliplr(p);

end

