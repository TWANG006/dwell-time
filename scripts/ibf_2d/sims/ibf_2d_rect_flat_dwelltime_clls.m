clc;
close all;
clear;

addpath(genpath('../../../ibf_engine'));

%% 0. Paramters
max_dt_diff = 50;    % max dwell time difference between consective P's
pixel_P_C_m = 2e-3;  

%% 1. Load the simulation result
data_dir = '../../../data/sims/';
surf_file = 'ibf_2d_70x30_0.12051mm_s50_std0.3.mat';
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
                 
interval_P = floor(pixel_P_C_m / pixel_m);

% X, Y coordinates of the dwell positions [m]
X_P_C = X(dw_range.y_s:interval_P:dw_range.y_e, dw_range.x_s:interval_P:dw_range.x_e);
Y_P_C = Y(dw_range.y_s:interval_P:dw_range.y_e, dw_range.x_s:interval_P:dw_range.x_e);


%% 4. CLLS Algorithm
% unused arguments
X_brf=0; 
Y_brf=0; 
Z_avg=0;

tic
[T_C, Z_removal_C, Z_residual_C, Z_to_remove_ca, Z_removal_ca_C, Z_residual_ca_C, C_C] = ...
    DwellTime2D_CLLS(Z_to_remove, X, Y, BRF_params, X_brf, Y_brf, Z_avg, X_P_C, Y_P_C, ca_range, max_dt_diff, 'model');
toc

% Fine sampling interval of dwell positions [pixel] 
pixel_P_F_m = 1e-3;                       
interval_P_1 = round(pixel_P_F_m / pixel_m);

% X, Y coordinates of the dwell positions [m]
X_P_F = X(dw_range.y_s:interval_P_1:dw_range.y_e, dw_range.x_s:interval_P_1:dw_range.x_e);
Y_P_F = Y(dw_range.y_s:interval_P_1:dw_range.y_e, dw_range.x_s:interval_P_1:dw_range.x_e);

% Interpolation
[T_F, Z_removal_F, Z_residual_F, Z_removal_ca_F, Z_residual_ca_F, C_F] = ...
    DwellTime2D_CLLS_Refine(Z_to_remove, X, Y, BRF_params, T_C, X_P_C, Y_P_C, X_P_F, Y_P_F);

%% 5. Show the results
DwellTime2D_Matrix_ShowResults(X, Y, Z_to_remove, Z_removal_C, Z_residual_C, Z_to_remove_ca, Z_removal_ca_C, Z_residual_ca_C, X_P_C, Y_P_C, T_C, ca_range);

%% 5. Save the coarse results
% coarse_filename = [surfaceMap_filename(1:end-4) '_CLLS_Coarse'];
% save([coarse_filename, '.mat'], ...
%     'BRF_params', ...
%     'x_start', 'y_start',...
%     'pixel_m', 'pixel_P_C_m',...
%     'X', 'Y', 'Z_to_remove', 'Z_removal_C', 'Z_residual_C',...
%     'C_C',...
%     'Z_to_remove_ca', 'Z_removal_ca_C', 'Z_residual_ca_C',...
%     'X_P_C', 'Y_P_C', 'T_C');
