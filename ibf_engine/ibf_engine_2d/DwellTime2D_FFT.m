function [T ... dwell map of the 2D IBF [s]
    , Z_removal... prediction of the removed height [m]
    , Z_residual ... prediction of the residual height [m]
    , T_dw... dwell map on the 2D dwell grid [s]
    , Z_to_remove_dw ...height to remove in the dwell grid [m]
    , Z_removal_dw...predection of the removal amount using T [m]
    , Z_residual_dw...prediction of the residual after removal using T [m]
    , Z_to_remove_ca...height to remove in the clear aperture [m] 
    , Z_removal_ca...predection of the removal amount using T in the clear aperture  [m]
    , Z_residual_ca...hprediction of the residual after removal using T in the clear aperture without detilt [m]
    ] = DwellTime2D_FFT ...  
    ( Z_to_remove... desired removal map
    , B... Calculated BRF
    , dw_range...dwell grid range
    , ca_range...clear aperture range
    )
% Purpose:
%   This function calculates the 2D dwell time map using the FFT method with
%   inverse filtering.
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

%% Calculate the dwell time by refining the inverse filtering threshold gamma
% 1. First time to get the gamma0
% Calculate T for the dwell grid
Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
Z_to_remove_dw = Z_to_remove_dw - nanmin(Z_to_remove_dw(:));
T_dw = DwellTime2D_FFT_InverseFilter(Z_to_remove_dw, B, 1, false);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_dw;

% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 

% Obtain the height to remove and height removal in the clear aperture
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
%Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

% Get gamma0
gamma0 = nanstd(Z_to_remove_ca(:), 1) / nanstd(Z_removal_ca(:), 1);
    
% Get optimized gamma
tic
gamma = DwellTime2D_FFT_OptimizeGamma(gamma0, Z_to_remove, Z_to_remove_dw, B, dw_range, ca_range, 'dwell', false);
toc

display([gamma0, gamma]);

% 2. Use the optimal gamma to do the computation again
T_dw = DwellTime2D_FFT_InverseFilter(Z_to_remove_dw, B, gamma, false);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_dw;

% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 

% Obtain the height to remove and height removal in the clear aperture
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    
% Get the Tdw on dwell grid
T_dw = T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

%% Obtain the entire aperture result
Z_residual = Z_to_remove - Z_removal;

%% Obtain the dwell grid result
Z_removal_dw = Z_removal(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
Z_residual_dw = Z_residual(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

%% Obtain the clear aperture results
Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca = Z_residual(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

% De-tilt
[x_ca,y_ca] = meshgrid(1:size(Z_to_remove_ca,2),1:size(Z_to_remove_ca, 1));
Z_to_remove_ca = RemoveSurface1(x_ca, y_ca, Z_to_remove_ca);
Z_to_remove_ca = Z_to_remove_ca - nanmin(Z_to_remove_ca(:));
Z_removal_ca = RemoveSurface1(x_ca, y_ca, Z_removal_ca);
Z_removal_ca = Z_removal_ca - nanmin(Z_removal_ca(:));
Z_residual_ca = RemoveSurface1(x_ca, y_ca, Z_residual_ca);

end