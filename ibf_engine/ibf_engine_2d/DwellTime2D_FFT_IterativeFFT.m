function [T, Z_removal, Z_residual...
    , T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca...
    ] = DwellTime2D_FFT_IterativeFFT...
    ( Z_to_remove... height to remove in the full aperture [m]
    , B...
    , dw_range, ca_range...
    , maxIters, PV_dif, RMS_dif, dwellTime_dif... stopping creteria
    )
% This function implement the Iterative DCT algorithm for dwell time
% calculation

%% Iteration 0
[T, Z_removal, Z_residual, T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, Z_residual_ca_woDetilt] =...
    DwellTime2D_FFT_IterativeFFT_OneIter(Z_to_remove, B, 0, dw_range, ca_range);

Z_residual_ca_pre = Z_residual_ca;
PV_pre = SmartPTV(Z_residual_ca(:));
% RMS_pre = nanstd(Z_residual_ca(:), 1);
dwellTime_pre = sum(T_dw(:));
num_iter = 1;

%% Iterations 1 to maxIters
while(true)
    % Calulcation of new dwell time
    [T, Z_removal, Z_residual, T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw, Z_to_remove_ca, Z_removal_ca, Z_residual_ca, Z_residual_ca_woDetilt] =...
        DwellTime2D_FFT_IterativeFFT_OneIter(Z_to_remove, B, - min(Z_residual_ca_woDetilt(:)), dw_range, ca_range);
            
    % Calculate stopping threshold
    PV_cur = SmartPTV(Z_residual_ca(:));
    %RMS_cur = nanstd(Z_residual_ca(:), 1);
    dwellTime_cur = sum(T_dw(:));
    
    if(num_iter>maxIters)
        display('Maximum number of iterations reached.');
        break;
%     elseif(abs(PV_cur-PV_pre )<PV_dif)
%         display('PV difference limit reached.');
%         break;
    elseif(nanstd(Z_residual_ca(:) - Z_residual_ca_pre(:), 1)<RMS_dif)
        display('RMS difference limit reached.');
        break;
    elseif(abs(dwellTime_pre-dwellTime_cur)<dwellTime_dif)
        display('Dwell time difference limit reached.');
        break;
    else
        Z_residual_ca_pre = Z_residual_ca;
        PV_pre = PV_cur;
        %RMS_pre = RMS_cur;
        dwellTime_pre = dwellTime_cur;
        num_iter = num_iter + 1;
    end    
end

end


function [T, Z_removal, Z_residual...
    , T_dw, Z_to_remove_dw, Z_removal_dw, Z_residual_dw...
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca, Z_residual_ca_woDetilt...
    ] = DwellTime2D_FFT_IterativeFFT_OneIter(Z_to_remove, B, hBias, dw_range, ca_range)

% The ca in dw range
ca_in_dw_y_s = ca_range.y_s - dw_range.y_s;
ca_in_dw_x_s = ca_range.x_s - dw_range.x_s;
ca_in_dw_y_e = ca_in_dw_y_s + ca_range.y_e - ca_range.y_s;
ca_in_dw_x_e = ca_in_dw_x_s + ca_range.x_e - ca_range.x_s;

Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
    
Z_to_remove = Z_to_remove - nanmin(Z_to_remove_ca(:));
Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) + hBias;

% 1. First time to get the gamma0
% Calculate T_dw for the dwell positions
T_dw = DwellTime2D_FFT_InverseFilter(Z_to_remove_dw, B, 1, false);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_dw;


% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 
Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_to_remove_ca = Z_to_remove_dw(ca_in_dw_y_s:ca_in_dw_y_e, ca_in_dw_x_s:ca_in_dw_x_e);
    
% Get gamma0
gamma0 = nanstd(Z_to_remove_ca(:), 1) / nanstd(Z_removal_ca(:), 1);
    
% Get optimized gamma
gamma = DwellTime2D_FFT_OptimizeGamma(gamma0, Z_to_remove, Z_to_remove_dw, B, dw_range, ca_range, 'dwell', false);
    
%  2. Use the optimal gamma to do the computation again
% Calculate T_dw for the dwell positions
T_dw = DwellTime2D_FFT_InverseFilter(Z_to_remove_dw, B, gamma, false);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_dw;

% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 
Z_to_remove_ca = Z_to_remove_dw(ca_in_dw_y_s:ca_in_dw_y_e, ca_in_dw_x_s:ca_in_dw_x_e);

%% Obtain the entire aperture result
% Calcualted Z_residual
Z_residual = Z_to_remove - Z_removal;

%% Obtain the dwell grid result
Z_removal_dw = Z_removal(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
Z_residual_dw = Z_residual(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

%% Obtain the clear aperture results
Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_residual_ca_woDetilt = Z_residual(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

% De-tilt
[x_ca,y_ca] = meshgrid(1:size(Z_to_remove_ca,2),1:size(Z_to_remove_ca, 1));
Z_to_remove_ca = RemoveSurface1(x_ca, y_ca, Z_to_remove_ca);
Z_to_remove_ca = Z_to_remove_ca - nanmin(Z_to_remove_ca(:));
Z_removal_ca = RemoveSurface1(x_ca, y_ca, Z_removal_ca);
Z_removal_ca = Z_removal_ca - nanmin(Z_removal_ca(:));
Z_residual_ca = RemoveSurface1(x_ca, y_ca, Z_residual_ca_woDetilt);

end