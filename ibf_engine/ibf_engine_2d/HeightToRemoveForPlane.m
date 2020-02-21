function Z_R = HeightToRemoveForPlane ...
    ( X ... x-coordinate [m]
    , Y ... y-coordinate [m]
    , Z ... height [m]
    )
% Function:
%   Z_R = HeightToRemoveForPlane(X, Y, Z)
%
% Purpose:
%   This function calculate the desired height to be removed for the plane.
%   The function first calculate the ROC R of the initial surface error,
%   only if the error is larger than 1e5 m will it be removed from the
%   desired height to be removed calculation.
%
% Inputs:
%   X: X meshgrid coordinates [m]
%   Y: Y meshgrid coordinates [m]
%   Z: Height error map [m]
%
% Outputs:
%   Z_R: Desired height to remove for the perfact plane
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

%% 1. Remove the 2D surface
[Z_R_ini, Z_f, f] = RemoveSurface2(X, Y, Z);

%% 2. Calculate the height to be removed
R = 1/f(5);

if(R > 1e5)
    Z_R = RemoveSurface1(X, Y ,Z);
    Z_R = Z_R - nanmin(Z_R(:));
else
    % Piston adjustment:
    Z = Z - nanmin(Z_R_ini(:));
    Z_R = Z - Z_f;
end

end