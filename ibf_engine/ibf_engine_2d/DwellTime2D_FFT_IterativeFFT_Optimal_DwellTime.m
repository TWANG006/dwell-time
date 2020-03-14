function [T, Z_removal, Z_residual...
    , dw_range, X_dw, Y_dw, T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca...
    ] = DwellTime2D_FFT_IterativeFFT_Optimal_DwellTime...
    ( Z_to_remove... height to remove in the full aperture [m]
    , X, Y...
    , B...
    , BRF_params...
    , ca_range...
    , maxIters, PV_dif, RMS_dif, dwellTime_dif... stopping creteria
    )

%% Get the BRF parameters
d_p = BRF_params.d_pix;         % diameter [pixel]
r_p = floor(0.5 * d_p);         % radius [pixel]

%% Compute the optical dwell time based on extention by radius of BRF 
% 1. Compute the optimal reference result based on r_p
dw_range.x_s = ca_range.x_s - r_p;   dw_range.x_e = ca_range.x_e + r_p;
dw_range.y_s = ca_range.y_s - r_p;   dw_range.y_e = ca_range.y_e + r_p;

[T, Z_removal, Z_residual...
    , T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca] = DwellTime2D_FFT_IterativeFFT(Z_to_remove, B, dw_range, ca_range, maxIters, PV_dif, RMS_dif, dwellTime_dif);

Opt_T_P_Sum = sum(T_dw(:));
Opt_Z_residual_ca = Z_residual_ca;
Z_residual_ca = Z_residual_ca*10;
% Opt_Z_residual_ca_RMS = nanstd(Z_residual_ca(:),1);
% Opt_Z_residual_ca_PV = SmartPTV(Z_residual_ca(:));

% 2. Iterate to find the better dwell time
s_p = floor(r_p/3);
Curr_T_P_Sum = Opt_T_P_Sum;
% Curr_Z_redisual_ca_RMS = Opt_Z_residual_ca_RMS + 2e-9;
% Curr_Z_residual_ca_PV = Opt_Z_residual_ca_PV + 2e-9;

while(s_p < r_p && ...
        (Curr_T_P_Sum >= Opt_T_P_Sum || ...
        nanstd(Z_residual_ca(:) - Opt_Z_residual_ca(:), 1) > 0.02e-9))
%         abs(Curr_Z_redisual_ca_RMS - Opt_Z_residual_ca_RMS) > 0.01e-9))... || ...
%         abs(Curr_Z_residual_ca_PV - Opt_Z_residual_ca_PV) > 0.01e-9))
    
    dw_range.x_s = ca_range.x_s - s_p;   dw_range.x_e = ca_range.x_e + s_p;
    dw_range.y_s = ca_range.y_s - s_p;   dw_range.y_e = ca_range.y_e + s_p;
    
    % Dwell grid coordinates
    X_dw = X(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    Y_dw = Y(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    
    [T, Z_removal, Z_residual...
    , T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca] = DwellTime2D_FFT_IterativeFFT(Z_to_remove, B, dw_range, ca_range, maxIters, PV_dif, RMS_dif, dwellTime_dif);
    
    s_p = s_p + 1;
    Curr_T_P_Sum = sum(T_dw(:));
%     Curr_Z_redisual_ca_RMS = nanstd(Z_residual_ca(:),1);
%     Curr_Z_residual_ca_PV = SmartPTV(Z_residual_ca(:));
end


end