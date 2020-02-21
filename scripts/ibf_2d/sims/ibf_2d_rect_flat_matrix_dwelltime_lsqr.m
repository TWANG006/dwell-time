clc;
close all;
clear;

addpath(genpath('../../../ibf_engine'));

%% 1. Load the simulation result
data_dir = '../../../data/sims/';
surf_file = 'ibf_2d_70x30_0.1mm_s50_std0.3.mat';
load([data_dir surf_file], 'X', 'Y', 'Z', 'Z_to_remove', 'pixel_m', 'zUnit_m');

%% 2. Define the BRF parameters
BRF_params.A = 0.4e-9; % Peak removal rate [m/s]
BRF_params.sigma_xy = FWHM2Sigma([2.7752e-3, 2.7752e-3]); % Sigma from FWHM
BRF_params.d = 10e-3; % Diameter
BRF_params.d_pix = round(BRF_params.d / pixel_m);
BRF_params.mu_xy = [0, 0]; % Center
r_p = floor(BRF_params.d_pix*0.5);

%% 3. Genreate the uniform machining path
% ca_range in [pixel% clear aperture (ca) range
x_range = 50e-3;                % range x [m]
y_range = 10e-3;                % range y [m]

x_s = BRF_params.d + pixel_m;   % x start [m]
y_s = BRF_params.d + pixel_m;   % x end [m]

x_e = x_s + x_range;            % x end [m]
y_e = y_s + y_range;            % y end [m]

ca_range.x_s = round(x_s / pixel_m); 
ca_range.x_e = round(x_e / pixel_m);
ca_range.y_s = round(y_s / pixel_m); 
ca_range.y_e = round(y_e / pixel_m);

% Get the dwell grid pixel range
dw_range.x_s = ca_range.x_s - r_p;   dw_range.x_e = ca_range.x_e + r_p;
dw_range.y_s = ca_range.y_s - r_p;   dw_range.y_e = ca_range.y_e + r_p;

pixel_P_m = 1e-3;                       
interval_P = floor(pixel_P_m / pixel_m);

% X, Y coordinates of the dwell positions [m]
X_P = X(dw_range.y_s:interval_P:dw_range.y_e, dw_range.x_s:interval_P:dw_range.x_e);
Y_P = Y(dw_range.y_s:interval_P:dw_range.y_e, dw_range.x_s:interval_P:dw_range.x_e);

%% 4. LSQR Algorithm
% unused arguments
X_brf=0; 
Y_brf=0; 
Z_avg=0;

tic;
[T, Z_removal, Z_residual, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, C] = ...
    DwellTime2D_LSQR(Z_to_remove, X, Y, BRF_params, X_brf, Y_brf, Z_avg, X_P, Y_P, ca_range, 0.35e-9, 'model');
toc;

% Z_test = Z_to_remove(y_start + 1:interval_P:end - y_start, x_start + 1:interval_P:end - x_start) + 1.56e-9;
% Z_test_max = max(Z_test(:));
% Si = (Z_test_max - Z_test) / Z_test_max;
% T = T.*Si;

%% 5. Show the LSQR results
DwellTime2D_Matrix_ShowResults(X, Y, Z_to_remove, Z_removal, Z_residual, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, X_P, Y_P, T, ca_range);