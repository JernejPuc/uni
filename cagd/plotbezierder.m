function db = plotbezierder (B, r, t, cfg)
% Opis :
%   plotbezier narice Bezierjevo krivuljo za dane kontrolne
%   tocke in seznam parametrov
%
% Definicija :
%   plotbezier (B,t)
%
% Vhodni podatki :
%   B    matrika velikosti n+1 x d, ki predstavlja kontrolne
%        tocke Bezierjeve krivulje stopnje n v
%        d- dimenzionalnem prostoru ,
%   t    seznam parametrov dolzine k, pri katerih racunamo
%        vrednost Bezierjeve krivulje

% Plt cfg
localHold = cfg(1);
showPoints = cfg(2);
showPoly = cfg(3);
showCurve = cfg(4);

% Local hold
if localHold
    hold on;
end

% Calc curve
db = bezierder(B, r, t);

% Get dims
deg = size(B, 1) - 1;
dim = size(db, 2);

% Calc ctrl pts
if r == 0
    cB = B;
else
    cB = diff(B, r) * prod(deg:-1:deg-r+1);
end

% Plt ctrl pts, poly and curve
if dim == 2
    if showPoints
        plot(cB(:, 1), cB(:, 2), 'o');
    end
    
    if showPoly
        for i = 1:size(cB, 1) - 1
            plot(cB(i:i + 1, 1), cB(i:i + 1, 2), 'k');
        end
    end
    
    if showCurve
        plot(db(:, 1), db(:, 2), 'k');
    end

elseif dim == 3
    if showPoints
        plot3(cB(:, 1), cB(:, 2), cB(:, 3), 'o');
    end
    
    if showPoly
        for i = 1:size(cB, 1) - 1
            plot3(cB(i:i + 1, 1), cB(i:i + 1, 2), cB(i:i + 1, 3), 'k');
        end
    end
    
    if showCurve
        plot3(db(:, 1), db(:, 2), db(:, 3));
    end
end

% Local hold
if localHold
    hold off;
end

end

