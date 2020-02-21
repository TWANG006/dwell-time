function [Z,Zf,f] = RemoveSurface2(X,Y,Z)
% Function:
%   [Z,Zf,f] = RemoveSurface2(X,Y,Z)
%
% Purpose:
%   This function removes the 2nd order curvature from the input surface
%   error map Z and returns the coefficients for those 2nd order variables.
%
% Inputs:
%   X: X meshgrid coordinates [m]
%   Y: Y meshgrid coordinates [m]
%   Z: Height error map [m]
%
% Outputs:
%   Z: Surface error map after removing curvature
%   Zf: Fitted surface error map with 2nd order polynomials
%   f: the fitting coefficients
%
% Info:
%   Contact: huanglei0114@gmail.com (Dr. HUANG Lei)
%   Copyright reserved.

idx = isfinite( Z(:) );
z = Z(idx);
x = X(idx);
y = Y(idx);

H = [ones(size(x)), ...
    x, y, ...
    x.^2, x.*y, y.^2,...
    ...x.^3, (x.^2).*y, x.*(y.^2), y.^3 ...
    ...x.^4, (x.^3).*y, (x.^2).*(y.^2), x.*(y.^3), y.^4 ...
    ];

f = pinv(H)*z;
% f = H\z;



Zf = (...
    f(1) ...
    + f(2)*X + f(3)*Y ...
    + f(4)*X.^2 + f(5)*X.*Y + f(6).*Y.^2 ...
    ...+ f(7)*X.^3 + f(8)*X.^2.*Y + f(9).*X.*Y.^2 + f(10).*Y.^3 ...
    ...+ f(11)*X.^4 + f(12)*X.^3.*Y + f(13)*X.^2.*Y.^2 + f(14).*X.*Y.^3 + f(15)*Y.^4 ...
    );

Z = Z - Zf;

end
