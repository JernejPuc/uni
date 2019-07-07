function rplotbezier(B, w, t, cfg)
% Plot rational bezier curve

% Plt cfg
localHold = cfg(1);
showPoints = cfg(2);
showFarin = cfg(3);
showPoly = cfg(4);
showCurve = cfg(5);

% Local hold
if localHold
    hold on;
end

% Calc curve
rb = rbezier(B, w, t);

% Get dims
dim = size(rb, 2);

% Plt ctrl pts, poly and curve
if dim == 2
    if showPoints
        plot(B(:, 1), B(:, 2), 'o');
    end

    if showFarin
        q = rfarin(B,w);
        plot(q(:, 1), q(:, 2), 'ro');
    end
    
    if showPoly
        for i = 1:size(B, 1) - 1
            plot(B(i:i + 1, 1), B(i:i + 1, 2), 'k');
        end
    end
    
    if showCurve
        plot(rb(:, 1), rb(:, 2));
    end

elseif dim == 3
    if showPoints
        plot3(B(:, 1), B(:, 2), B(:, 3), 'o');
        
    end

    if showFarin
        q = rfarin(B,w);
        plot3(q(:, 1), q(:, 2), q(:, 3), 'ro');
    end
    
    if showPoly
        for i = 1:size(B, 1) - 1
            plot3(B(i:i + 1, 1), B(i:i + 1, 2), B(i:i + 1, 3), 'k');
        end
    end
    
    if showCurve
        plot3(rb(:, 1), rb(:, 2), rb(:, 3));
    end
end

% Local hold
if localHold
    hold off;
end

end

