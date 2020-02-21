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
BRF_params.lat_res_brf = pixel_m;
d_pix = BRF_params.d_pix;
r_pix = d_pix * 0.5;

options = struct(...
    'maxIters', 10000,...
    'RMS_dif', 0.04e-9,...[m]
    'dwellTime_dif', 1e-2,...[s]
    'lambda', 2e-3,... % parameter for the total variation (TV) penalty
    'isDownSampling', true,...
    'samplingInterval', 1e-3...
    );
%     'alpha', 1,... % step size for the additive iterations6
  

%% 3. Genreate the uniform machining path
% ca_range in [pixel% clear aperture (ca) range
x_range = 50e-3;                % range x [m]
y_range = 5e-3;                % range y [m]

x_s = BRF_params.d + pixel_m;   % x start [m]
y_s = BRF_params.d + 2.5e-3 + pixel_m;   % x end [m]

x_e = x_s + x_range;            % x end [m]
y_e = y_s + y_range;            % y end [m]

ca_range.x_s = round(x_s / pixel_m); 
ca_range.x_e = round(x_e / pixel_m);
ca_range.y_s = round(y_s / pixel_m); 
ca_range.y_e = round(y_e / pixel_m);

%% 4. Algorithm
% unused arguments
X_brf=0; 
Y_brf=0; 
Z_avg=0;

tic;
[ B, X_B, Y_B...
, X, Y, Z_removal, Z_residual...
, T_P, X_P, Y_P...
, X_dw, Y_dw...
, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...dwell grid results [m]
, X_ca, Y_ca, Z_to_remove_ca, Z_removal_ca, Z_residual_ca...
    ] = DwellTime2D_Bayesian...
    ( Z_to_remove... height to remove [m]
    , 'model'...
    , BRF_params... BRF parameters
    , X_brf, Y_brf, Z_avg...
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , options...
    );
toc;

%% Display results
DwellTime2D_FFT_ShowResults(...
    B, X_B, Y_B,...
    X, Y,  Z_to_remove, Z_removal, Z_residual,... 
    T_P,  X_P, Y_P,...
    X_dw, Y_dw,  Z_to_remove_dw, Z_removal_dw, Z_residual_dw,...
    X_ca, Y_ca, Z_to_remove_ca, Z_removal_ca, Z_residual_ca);