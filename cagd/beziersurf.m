function beziersurf(Bx, By, Bz, u, v, cfg)
%BEZIERSURF Summary of this function goes here
%   Detailed explanation goes here

% Optional config
if nargin < 6
    localHold = true;
    showGrid = true;
    showSurf = true;
else
    localHold = cfg(1);
    showGrid = cfg(2);
    showSurf = cfg(3);
end

if localHold
    hold on;
end

% Ctrl poly grid + mesh settings
if showGrid
    m = mesh(Bx, By, Bz);
    set(m, 'EdgeColor', 'k');
    set(m, 'FaceAlpha', 0);
    set(m, 'LineWidth', 1);
end

% Calc + surf
if showSurf
    [bx, by, bz] = bezier2(Bx, By, Bz, u, v);
    surf(bx, by, bz);
end

grid on;

if localHold
    hold off;
end

end

