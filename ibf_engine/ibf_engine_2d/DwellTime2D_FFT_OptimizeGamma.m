function gamma = DwellTime2D_FFT_OptimizeGamma(gamma0, Z_to_remove, Z_to_remove_dw, B, dw_range, ca_range, flag, useDCT)
% Optimize the gamma parameter used in the inverse filtering algorithm by
% minimizing: f(gamma) = RMS(Z_residual_ca)
% Inputs:
%           gamma0: initial guess of the gamma
%      Z_to_remove: the height to remove in the entire aperture [m]
%   Z_to_remove_dw: the height to remove in the dwell gird [m]
%                B: the BRF of the ion gun
%                d: the diameter of B
%                r: the radius of b
%             flag: 'entire' for entire aperture calculation; 'dwell' for dwell grid calculation

options = optimoptions('patternsearch','Display', 'off');

if(strcmp(flag, 'entire'))
    fun = @(gamma)Objective_Func_Entire(gamma, Z_to_remove, B, dw_range, ca_range, useDCT);
    gamma = patternsearch(fun, gamma0, [], [], [], [], [], [], [], options);
elseif(strcmp(flag, 'dwell'))
    fun = @(gamma)Objective_Func_DwellGrid(gamma, Z_to_remove, Z_to_remove_dw, B, dw_range, ca_range, useDCT);
    gamma = patternsearch(fun, gamma0, [], [], [], [], [], [], [], options); 
else
    gamma = gamma0;
end

end

%% Merit functions fo f(gamma)
function fGamma = Objective_Func_Entire(gamma, Z_to_remove, B, dw_range, ca_range, useDCT)
% Calculate T_dw for the dwell positions
Tms = DwellTime2D_FFT_InverseFilter(Z_to_remove, B, gamma, useDCT);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = Tms(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
            
% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 

% Calculate the residual
Z_residual_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e) - Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

fGamma = nanstd(Z_residual_ca(:), 1);

end

function fGamma = Objective_Func_DwellGrid(gamma, Z_to_remove, Z_to_remove_dw, B, dw_range, ca_range, useDCT)

% The ca in dw range
ca_in_dw_y_s = ca_range.y_s - dw_range.y_s + 1;
ca_in_dw_x_s = ca_range.x_s - dw_range.x_s + 1;
ca_in_dw_y_e = ca_in_dw_y_s + ca_range.y_e - ca_range.y_s;
ca_in_dw_x_e = ca_in_dw_x_s + ca_range.x_e - ca_range.x_s;

% Calculate T_dw for the dwell positions
T_dw = DwellTime2D_FFT_InverseFilter(Z_to_remove_dw, B, gamma, useDCT);

% Only keep T in the dwell grid and let others to be 0
T = zeros(size(Z_to_remove));
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_dw;
            
% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 

% Calculate the residual
Z_residual_ca = Z_to_remove_dw(ca_in_dw_y_s:ca_in_dw_y_e, ca_in_dw_x_s:ca_in_dw_x_e)...
    - Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

fGamma = nanstd(Z_residual_ca(:), 1);

end