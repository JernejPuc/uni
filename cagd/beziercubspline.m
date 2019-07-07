function B = beziercubspline(u,D)
% Opis:
%   beziercubspline izracuna sestavljeno Bezierjevo krivuljo
%   stopnje 3, ki je dvakrat zvezno odvedljiva v stikih
%
% Definicija:
%   B = beziercubspline(u,D)
%
% Vhodna  podatka:
%   u   seznam parametrov delitve dolzine m+1,
%   D   matrika, v kateri vsaka izmed m+3 vrstic predstavlja
%       eno kontrolno tocko sestavljene  krivulje
%
% Izhodni  podatek:
%   B   seznam dolcine m, v kateri je vsak element matrika s
%       stirimi vrsticami, ki dolocajo kontrolne tocke kosa
%       sestavljene krivulje

% Get dims
mp1 = size(u, 2);
d = size(D, 2);

% Get m
m = mp1 - 1;

% Allocate
B = cell(1, m);

for i = 1:m
    B{i} = zeros(4, d);
end

% Diference
Del = diff(u);

% Formule na robu (b_{0,1,2}^{1}, b_{3,2,1}^{m})
B{1}(1, :) = D(1, :);
B{1}(2, :) = D(2, :);
B{1}(3, :) = (Del(2) / (Del(1) + Del(2))) * D(2, :);
B{1}(3, :) = B{1}(3, :) + (Del(1) / (Del(1) + Del(2))) * D(3, :);

B{m}(4, :) = D(m+3, :);
B{m}(3, :) = D(m+2, :);
B{m}(2, :) = (Del(m) / (Del(m-1) + Del(m))) * D(m+1, :);
B{m}(2, :) = B{m}(2, :) + (Del(m-1) / (Del(m-1) + Del(m))) * D(m+2,:);

% Za i = 1, ..., m-2
for i = 2:m-1
    B{i}(2, :) = ((Del(i) + Del(i+1)) / (Del(i-1) + Del(i) + Del(i+1))) * D(i+1, :);
    B{i}(2, :) = B{i}(2, :) + (Del(i-1) / (Del(i-1) + Del(i) + Del(i+1))) * D(i+2, :);
    
    B{i}(3, :) = (Del(i+1) / (Del(i-1) + Del(i) + Del(i+1))) * D(i+1, :);
    B{i}(3, :) = B{i}(3, :) + ((Del(i) + Del(i-1)) / (Del(i-1) + Del(i) + Del(i+1))) * D(i+2, :);
end

% Za i = 1, ..., m-1
for i = 2:m
    B{i-1}(4, :) = Del(i) / (Del(i) + Del(i-1)) * B{i-1}(3, :);
    B{i-1}(4, :) = B{i-1}(4, :) + Del(i-1) / (Del(i) + Del(i-1)) * B{i}(2, :);
    B{i}(1, :) = B{i-1}(4, :);
end

end

