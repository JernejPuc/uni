function R = matrotate(phi, mode)
%AFFINETRANS Summary of this function goes here
%   Detailed explanation goes here

% To radians
if strcmp(mode, 'torad')
    phi = phi * pi/180;
end

% Construct rotation matrix
if size(phi, 2) == 1
    % Dim = 2
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
else
    % Dim = 3
    phi_x = phi(1);
    phi_y = phi(2);
    phi_z = phi(3);
    
    Rx = [1 0 0; 0 cos(phi_x) -sin(phi_x); 0 sin(phi_x) cos(phi_x)];
    Ry = [cos(phi_y) 0 sin(phi_y); 0 1 0; -sin(phi_y) 0 cos(phi_y)];
    Rz = [cos(phi_z) -sin(phi_z) 0; sin(phi_z) cos(phi_z) 0; 0 0 1];
    
    R = Rz * Ry * Rx;
end

end

