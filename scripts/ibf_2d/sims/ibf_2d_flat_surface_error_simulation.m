clc;
close all;
clear;

addpath(genpath('../../../ibf_engine'));

%% 0. Define the pixel size and height units
pixel_m = 0.1e-3;   % resolution [m]
zUnit_m = 50e-9;    % height [m]
std_m = 0.3e-9;     % noise level

%% 1. Simulate the Surface error map using Legendre polynomials

range_x_mm = 100;  % mm
range_y_mm = 30;   % mm

range_x = range_x_mm*1e-3;   
range_y = range_y_mm*1e-3;   

x_mm = 0:pixel_m:range_x-pixel_m;
y_mm = 0:pixel_m:range_y-pixel_m;

current_mfile = ['../../../data/sims/ibf_2d_' num2str(range_x_mm) 'x' num2str(range_y_mm) '_' num2str(pixel_m*1e3) 'mm' '_s' num2str(zUnit_m*1e9) '_std' num2str(std_m*1e9)];

% Define the number of pixels along x and y [pixel]
x = round(x_mm / pixel_m);      % x range [pixel]
y = round(y_mm / pixel_m);      % y range [pixel]

% normalize x, y to [-1, 1]
xx = -1 + 2.*(x - min(x))./(max(x) - min(x)); 
yy = -1 + 2.*(y - min(y))./(max(y) - min(y));

% Use Lengendre polynomials to simulate the surface map
Z = LegendreP2D_XYJC(xx, yy, [4,5,6,7,8,9,10], [-1,0,-1,2,0,-1, -0.5]) * zUnit_m;

% Add noise
GN = randn(size(Z))*std_m;  % Gaussian noise with 0 mean and 4 variance
Z = Z + GN;

% Generate X, Y coordinate meshgrid
[X, Y] = meshgrid(x, y);
X = X*pixel_m;
Y = Y*pixel_m;

% Show the simulated surface map
ShowSurfaceMap(X, Y, Z, 'Simulated Surface Profile');


%% 2. Calculate the Height to be removed Z_R
% Remove the 2D tilt and make Z nonnegative
% Z = remove_2D_tilt(X, Y, Z);
% Z = Z - nanmin(Z(:));

% Calculate and show the height to be removed
Z_to_remove = HeightToRemoveForPlane(X, Y, Z);
ShowSurfaceMap(X, Y, Z_to_remove, 'Height to remove');

%% 3. Save the simultaed surface map and the height to remove
save([current_mfile, '.mat'], 'X', 'Y', 'Z', 'Z_to_remove', 'pixel_m', 'zUnit_m');

