clc;
close all;
clear;

addpath(genpath('../../../ibf_engine'));

%% 1. Load the simulation result
data_dir = '../../../data/sims/';
surfaceMap_filename = 'ibf_2d_70x30_0.12051mm_s50_std0.3.mat';
load(surfaceMap_filename, 'X', 'Y', 'Z', 'Z_to_remove', 'pixel_m', 'zUnit_m');

%% 2. Define the BRF parameters
BRF_params.A = 0.4e-9; % Peak removal rate [m/s]
BRF_params.sigma_xy = FWHM2Sigma([2.7752e-3, 2.7752e-3]); % Sigma from FWHM
BRF_params.d = 10e-3; % Diameter
BRF_params.d_pix = round(BRF_params.d / pixel_m);
BRF_params.mu_xy = [0, 0]; % Center

%% 3. Genreate the uniform machining path
% Number of boundary pixels along x and y [pixel]
x_start = ceil(BRF_params.d * 0.5 / pixel_m);
y_start = ceil(BRF_params.d * 0.5 / pixel_m);

pixel_P_m = 1e-3;                       
interval_P = floor(pixel_P_m / pixel_m);

% X, Y coordinates of the dwell positions [m]
X_P = X(y_start + 1:interval_P:end - y_start, x_start + 1:interval_P:end - x_start);
Y_P = Y(y_start + 1:interval_P:end - y_start, x_start + 1:interval_P:end - x_start);

%% 4. TSVD Algorithm
tic;
[T, Z_removal, Z_residual, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, C] = ...
    DwellTime2D_TSVD(Z_to_remove, X, Y, BRF_params, X_P, Y_P, 0.33e-9);
toc;

%% 5. Show the TSVD results
DwellTime2D_Matrix_ShowResults(X, Y, Z_to_remove, Z_removal, Z_residual, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, X_P, Y_P, T, BRF_params, pixel_P_m);