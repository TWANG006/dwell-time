clc;
close all;
clear;

addpath(genpath('../../../ibf_engine'));

%% 1. Load the simulation result
data_dir = '../../../data/sims/';
surf_file = 'ibf_2d_70x30_0.12051mm_s50_std0.3.mat';
load([data_dir surf_file], 'X', 'Y', 'Z', 'Z_to_remove', 'pixel_m', 'zUnit_m');

%% 2. Define the BRF parameters
BRF_params.A = 1e-9;                                      % Peak removal rate [m/s]
BRF_params.sigma_xy = FWHM2Sigma([2.7752e-3, 2.7752e-3]); % Sigma from FWHM
BRF_params.d_pix = 84;                                    % Diameter [pix]  
BRF_params.d = BRF_params.d_pix * pixel_m;                % Diameter [m]
BRF_params.mu_xy = [0, 0];                                % Center [m]
BRF_params.lat_res_brf = pixel_m;
d_pix = BRF_params.d_pix;
r_pix = d_pix * 0.5;

% 2.0 Inverse filtering without iterative updates
options = struct(...
    'Algorithm', 'FFT',...
    'maxIters', 10, ...
    'PV_dif', 0.001e-9, ...[m]
    'RMS_dif', 0.02e-9, ...[m]
    'dwellTime_dif', 60, ...[s]
    'isDownSampling', false, ...
    'samplingInterval', 1e-3 ... [m]
);

% 2.1 Use Inner iteration only
% options = struct(...
%     'Algorithm', 'Iterative-FFT',...
%     'maxIters', 10, ...       
%     'PV_dif', 0.001e-9, ...[m]
%     'RMS_dif', 0.02e-9, ...[m]
%     'dwellTime_dif', 60, ...[s]
%     'isDownSampling', false, ...  wheter downsample or not
%     'samplingInterval', 1e-3 ... [m]
% );

% 2.2 Use both inner and outer iterations
% options = struct(...
%     'Algorithm', 'Iterative-FFT-Optimal-DwellTime',...
%     'maxIters', 10, ...
%     'PV_dif', 0.001e-9, ...[m]
%     'RMS_dif', 0.02e-9, ...[m]
%     'dwellTime_dif', 60, ...[s]
%     'isDownSampling', true, ...
%     'samplingInterval', 1e-3 ... [m]
% );

%% 3. Define the clear aperture 

X_brf = 0;
Y_brf = 0;
Z_avg = 0;
brf_mode = 'model';
tmin = 0;
tmax = 1;

[m, n] = size(Z_to_remove);
ca_range.x_s = BRF_params.d_pix + 1; ca_range.x_e = n - BRF_params.d_pix;
ca_range.y_s = BRF_params.d_pix + 1; ca_range.y_e = m - BRF_params.d_pix;


%% 4. Run the algorithm
tic;
[ B, X_B, Y_B...
, X, Y, Z_removal, Z_residual...
, T_P, T_P_Real, X_P, Y_P...
, X_dw, Y_dw, dw_range...
, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...dwell grid results [m]
, X_ca, Y_ca, Z_to_remove_ca, Z_removal_ca, Z_residual_ca...
    ] = DwellTime2D_FFT_Full...
    ( Z_to_remove... height to remove [m]
    , BRF_params... BRF parameters
    , brf_mode...
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin, tmax...
    , options...
    );
toc;

%% 5. Display results
DwellTime2D_FFT_ShowResults(...
    B, X_B, Y_B,...
    X, Y,  Z_to_remove, Z_removal, Z_residual,... 
    T_P,  X_P, Y_P,...
    X_dw, Y_dw,  Z_to_remove_dw, Z_removal_dw, Z_residual_dw,...
    X_ca, Y_ca, Z_to_remove_ca, Z_removal_ca, Z_residual_ca);