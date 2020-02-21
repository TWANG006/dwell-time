function [T_F ... refined dwell map of the 2D IBF [s]
    , z_removal_F... prediction of the removed height
    , z_residual_F ... prediction of the residual height
    , Z_removal_ca_F...predection of the removal amount using T
    , Z_residual_ca_F...prediction of the residual after removal using T
    , C_F... the matrix C for the refinement
    ] = DwellTime2D_CLLS_Refine(Z_to_remove... height to be removed [m]
    , X... X grid points of the z_to_remove
    , Y... Y grid points of the z_to_remove
    , BRF_params... BRF parameters
    , X_brf, Y_brf, Z_avg...
    , T... oritingal dwell time [s]
    , X_P... x coordinates of the dwell points, size = m x n
    , Y_P... y coordinates of the dwell points, size = m x n
    , X_P_F... 
    , Y_P_F...
    , ca_range...
    , resampling_method...
    )
% Purpose:
%   This function resample the dwell tiem map calculated using CLLS to the
%   desired sampling interval
%
% Reference:
%   Wang, T., Huang, L., Vescovi, M., Kuhne, D., Tayabaly, K., 
%   Bouet, N., & Idir, M. (2019). Study on an effective one-dimensional 
%   ion-beam figuring method. Optics express, 27(11), 15368-15381.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.
%% Build the dwell position list
P_F = [X_P_F(:), Y_P_F(:)];

Nt_F = size(P_F, 1);
Nr = numel(Z_to_remove);
scaleFactor = (X_P_F(1, 2) - X_P_F(1, 1)) / (X_P(1, 2) - X_P(1,1));
scaleFactor = scaleFactor.^2;
d_p = BRF_params.d_pix;

%% Build the refined matrix C_F for the calculation
[C_F, d_F, C_T_F] = DwellTime2D_Assemble_C_d(Nr, Nt_F, BRF_params, Z_to_remove, X, Y, P_F, X_brf, Y_brf, Z_avg, ca_range, resampling_method);

%% Interpolation to construct finer T_F from T
T_F = interp2(X_P, Y_P, T, X_P_F, Y_P_F, 'spline') * scaleFactor;
T_F = T_F(:);

% Clear aperture results
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_removal_ca_F = C_F * T_F;
Z_residual_ca_F = d_F - Z_removal_ca_F;

Z_removal_ca_F = reshape(Z_removal_ca_F, size(Z_to_remove_ca));
Z_removal_ca_F = RemoveSurface1(X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Z_removal_ca_F);
Z_removal_ca_F = Z_removal_ca_F - nanmin(Z_removal_ca_F(:));

Z_residual_ca_F = reshape(Z_residual_ca_F, size(Z_to_remove_ca));
Z_residual_ca_F = RemoveSurface1(X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Z_residual_ca_F);



% Full results
z_removal_F = C_T_F * T_F;
z_residual_F = Z_to_remove(:) - z_removal_F;
z_removal_F = reshape(z_removal_F, size(Z_to_remove));
z_residual_F = reshape(z_residual_F, size(Z_to_remove));

T_F = reshape(T_F, size(X_P_F));

end