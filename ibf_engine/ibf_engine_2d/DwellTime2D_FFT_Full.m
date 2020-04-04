function [B... BRF
    , X_B, Y_B...BRF Coordinates
    , X, Y...full aperture coordinates [m]
    , Z_removal, Z_residual... full aperture results [m]
    , T_P... dwell time on the dwell grid [s]
    , T_P_Real...
    , X_P, Y_P...dwell grid 
    , X_dw, Y_dw, dw_range ... dwell grid coordinates [m]
    , Z_to_remove_dw, Z_removal_dw, Z_residual_dw...dwell grid results [m]
    , X_ca, Y_ca... clear aperture coordinates [m]
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca...[m]
    ] = DwellTime2D_FFT_Full...
    ( Z_to_remove... height to remove [m]
    , BRF_params... BRF parameters
    , BRF_mode... 'avg' or 'model'
    , X_BRF, Y_BRF, Z_BRF ... averaged z & its coords
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , tmin, tmax...
    , options...
    , ratio...
    )
% Purpose:
%   This function implement the optimized FFT-based dwell time calculation
%   (deconvolution) algorithm with inverse filtering and the RIFTA. 
%
% Reference:
%   [1] T. Wang, L. Huang, H. Kang, H. Choi, D. W. Kim, K. Tayabaly,
%   and M. Idir, “Rifta: a robust iterative fourier transform-based 
%   dwell time algo-rithm for ion beam figuring,” 
%   Sci. Reports, Under peer review(2019)
%
%   [2] Wang, T., Huang, L., Tayabaly, K., & Idir, M. (2019, November). 
%   Study on the performances of dwell time algorithms in ion beam figuring. 
%   In Optifab 2019 (Vol. 11175, p. 111750M). International Society for 
%   Optics and Photonics.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

%% 0. Set the default options for the function
defaultOptions = struct(...
    'Algorithm', 'Iterative-FFT',...
    'maxIters', 10,...
    'PV_dif', 0.001e-9,...[m]
    'RMS_dif', 0.001e-9,...[m]
    'dwellTime_dif', 60,...[s]
    'isDownSampling', false...
);

%% 1. Deal with input arguments
% check invalid method
if nargin == 10
    options = defaultOptions; 
end
if nargin == 11
    ratio = 1; 
end
 
%% 2. Construct the BRF using the BRF parameters
% Release the BRF parameters
A = BRF_params.A;               % peak removal rate [m/s]
sigma_xy = BRF_params.sigma_xy; % standard deviation [m]
mu_xy = [0,0];                  % center is 0 [m]
brf_res = BRF_params.lat_res_brf;
brf_pix = BRF_params.d_pix;

brf_d = brf_res*brf_pix;    % [m] range of Measurement_Area
brf_r = brf_d * 0.5;        % radius of the brf

[X_B, Y_B] = meshgrid(-brf_r:pixel_m:brf_r, -brf_r:pixel_m:brf_r);


% X, Y coordinates for B
% ratio_B = pixel_m / BRF_params.lat_res_brf;
%     
% % Resample the BRF grid
% X_B =  imresize(X_BRF, 1/ratio_B);
% Y_B =  imresize(Y_BRF, 1/ratio_B);

% 
% [X_B, Y_B] = meshgrid(-r_p : r_p - 1, -r_p : r_p - 1);        
% X_B = X_B * pixel_m;
% Y_B = Y_B * pixel_m;

% Get B
if strcmpi(BRF_mode, 'avg')
    B = interp2(X_BRF, Y_BRF, Z_BRF, X_B, Y_B, 'spline');
else
    B = BRFGaussian2D(X_B, Y_B, 1, [A, sigma_xy, mu_xy]); 
end

d_p = size(B,1);         % diameter [pixel]
r_p = floor(0.5 * d_p);  % radius [pixel]

% reset the BRF params
BRF_params.lat_res_brf = pixel_m;
BRF_params.d_pix = d_p;
    
%% 3. Define the dwell grid
% Get the size of the full aperture
[mM, nM] = size(Z_to_remove);   

% Get the dwell grid pixel range
dw_range.x_s = ca_range.x_s - r_p - 1;   dw_range.x_e = ca_range.x_e + r_p + 1;
dw_range.y_s = ca_range.y_s - r_p - 1;   dw_range.y_e = ca_range.y_e + r_p + 1;

% Determine if the range is valid
if(dw_range.x_s<1 || dw_range.x_e >nM || dw_range.y_s <1 || dw_range.y_e > mM)
    error(['Invalid clear aperture range with [' num2str(dw_range.x_s) ', ' num2str(dw_range.x_e) ']' ' and ' '[' num2str(dw_range.y_s) ', ' num2str(dw_range.y_e) ']']);
else
    % Get X, Y pixel coordinates
    [X, Y] = meshgrid(0:nM -1, 0:mM - 1);
    
    % Transfer to physical coordinates
    X = X * pixel_m;    
    Y = Y * pixel_m;
    
    % Dwell grid coordinates
    X_dw = X(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    Y_dw = Y(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    
    % Clear aperture coordinates
    X_ca = X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
    Y_ca = Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);    
end

%% 4. Real FFT algorithm
% Iterative FFT on dwell grid
if strcmp(options.Algorithm, 'Iterative-FFT')
    maxIters = options.maxIters;
    PV_dif = options.PV_dif;
    RMS_dif = options.RMS_dif;
    dwellTime_dif = options.dwellTime_dif;
    
    [~, Z_removal, Z_residual...
    , T_P, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca] = DwellTime2D_FFT_IterativeFFT(Z_to_remove, B, dw_range, ca_range, maxIters, PV_dif, RMS_dif, dwellTime_dif);
elseif strcmp(options.Algorithm, 'Iterative-FFT-Optimal-DwellTime')
    maxIters = options.maxIters;
    PV_dif = options.PV_dif;
    RMS_dif = options.RMS_dif;
    dwellTime_dif = options.dwellTime_dif; 
    
    [~, Z_removal, Z_residual...
     , dw_range, X_dw, Y_dw, T_P, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
     , Z_to_remove_ca, Z_removal_ca, Z_residual_ca] = DwellTime2D_FFT_IterativeFFT_Optimal_DwellTime(Z_to_remove, X, Y, B, BRF_params, ca_range, maxIters, PV_dif, RMS_dif, dwellTime_dif);
elseif strcmp(options.Algorithm, 'FFT')
    [~, Z_removal, Z_residual...
    , T_P, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca] = DwellTime2D_FFT(Z_to_remove, B, dw_range, ca_range);
else
    error('Invalid FFT algorithm chosen. Should be either Iterative-FFT or FFT');
end

T_P = T_P * ratio;

%% 5. Downsampling if used
% Use a sparser dwell grid
if options.isDownSampling == true
    % Obtain the sampling interval
    pixel_P_m = options.samplingInterval;
    interval_P_m = pixel_P_m / pixel_m;
    
    % Down sample the dwell grid
    X_P =  imresize(X_dw, 1/interval_P_m);
    Y_P =  imresize(Y_dw, 1/interval_P_m);
    
    % Dump X_P, Y_P & X_P_F, Y_P_F dwell point positions into a 2D array as
    %   |  u1    v1 |   P1
    %   |  u2    v2 |   P2
    %   | ...   ... |  ...
    %   |  uNt   vNt|   PNt
    P = [X_P(:), Y_P(:)];
    
    % Get the numbers of IBF machining points and sampling points of the surface error map R
    Nt = size(P, 1);
    Nr = numel(Z_to_remove);
    
    % Assemble the BRF matrix C, size(C) = Nr x Nt and vector d
    [C, d, C_T] = DwellTime2D_Assemble_C_d(Nr, Nt, BRF_params, Z_to_remove, X, Y, P, X_BRF, Y_BRF, Z_BRF, ca_range, BRF_mode);

    % Down sample T_dw
    T_P = imresize(T_P, 1/interval_P_m, 'bicubic') * interval_P_m.^2;
    T_P_Real = T_P + tmin*ceil(max(T_P(:)) / tmax);
    T_P_v = T_P_Real(:);

    % Clear aperture results
    Z_removal_ca = C * T_P_v;
    Z_residual_ca = d - Z_removal_ca;
    Z_removal_ca = reshape(Z_removal_ca, size(Z_to_remove_ca));
    Z_residual_ca = reshape(Z_residual_ca, size(Z_to_remove_ca));
    Z_to_remove_ca = reshape(d, size(Z_to_remove_ca));
    
    % Detilt
    Z_to_remove_ca = RemoveSurface1(X_ca, Y_ca, Z_to_remove_ca);
    Z_to_remove_ca = Z_to_remove_ca - nanmin(Z_to_remove_ca(:));
    Z_removal_ca = RemoveSurface1(X_ca, Y_ca, Z_removal_ca);
    Z_removal_ca = Z_removal_ca - nanmin(Z_removal_ca(:));
    Z_residual_ca = RemoveSurface1(X_ca, Y_ca, Z_residual_ca);
    
    % Full aperture results
    Z_removal = C_T * T_P_v;
    Z_residual = Z_to_remove(:) - Z_removal;
    Z_residual = Z_residual - nanmean(Z_residual);
    Z_removal = reshape(Z_removal, size(Z_to_remove));
    Z_residual = reshape(Z_residual, size(Z_to_remove));  
    
    % Dwell grid results
    Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    Z_removal_dw = Z_removal(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    Z_residual_dw = Z_residual(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);   
else
    X_P = X_dw;
    Y_P = Y_dw;   
    T_P_Real = T_P;
end

end