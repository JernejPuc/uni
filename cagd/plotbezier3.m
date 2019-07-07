function plotbezier3(Bx,By,Bz,k,cfg)
%PLOTBEZIER3 Summary of this function goes here
%   Detailed explanation goes here

% Optional config
if nargin < 4
    k = 10;
end

if nargin < 5
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
    [TRI,U] = trimeshgrid(k);
    b = bezier3(Bx,By,Bz,U);
    trisurf(TRI,b(:,1),b(:,2),b(:,3))
end

grid on;

if localHold
    hold off;
end

end

