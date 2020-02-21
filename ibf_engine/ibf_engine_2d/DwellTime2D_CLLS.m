function [T ... [s] dwell map of the 2D IBF
    , z_removal... prediction of the removed height
    , z_residual ... prediction of the residual height
    , Z_to_remove_ca...
    , z_removal_ca...predection of the removal amount using T
    , z_residual_ca...prediction of the residual after removal using T
    , C ... matrix C resembled from BRF
    ] = DwellTime2D_CLLS ...  
    ( Z_to_remove... desired removal map
    , X... X grid points of the z_to_remove
    , Y... Y grid points of the z_to_remove
    , BRF_params... BRF parameters
    , X_brf, Y_brf, Z_avg...
    , X_P... x coordinates of the dwell points, size = m x n
    , Y_P... y coordinates of the dwell points, size = m x n
    , ca_range...
    , max_dt_diff... max dwell time diff constraint [s]     
    , resampling_method...
    )
% Purpose:
%   This function implements the matrix-based CLLS dwell time algorithm
%
%                 Nt
%   R(xk, yk) = Sigma C(xk - ui, yk-vi)T(ui, vi)
%                 i=1
%
% Reference:
%   Wang, T., Huang, L., Vescovi, M., Kuhne, D., Tayabaly, K., 
%   Bouet, N., & Idir, M. (2019). Study on an effective one-dimensional 
%   ion-beam figuring method. Optics express, 27(11), 15368-15381.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

%% Dump X_P, Y_P & X_P_F, Y_P_F dwell point positions into a 2D array as
%   |  u1    v1 |   P1
%   |  u2    v2 |   P2
%   | ...   ... |  ...
%   |  uNt   vNt|   PNt
P = [X_P(:), Y_P(:)];

%% Get the numbers of IBF machining points and sampling points of the surface error map R
Nt = size(P, 1);
Nr = numel(Z_to_remove);

%% Assemble the BRF matrix C, size(C) = Nr x Nt and vector d
[C, d, C_T] = DwellTime2D_Assemble_C_d(Nr, Nt, BRF_params, Z_to_remove, X, Y, P, X_brf, Y_brf, Z_avg, ca_range, resampling_method);

%% Assemble the constraints matrix A and vector b
[m, n] = size(X_P);
[A, b] = DwellTime2D_CLLS_Assemble_A_b(m, n, max_dt_diff);
% [A, b] = DwellTime2D_CLLS_Assemble_A_b_smart(m, n, Z_to_remove, 0.002);

%% Assemble LB and UB for t
LB = zeros(Nt, 1);
UB = Inf(Nt, 1);

%% CLLS algorithm to calculate the coarse T
OptOptions = optimset( 'MaxIter', 5e6, 'Algorithm','active-set');
% T = lsqlin(C, d, A, b, [], [], LB, UB, [], OptOptions);
T = lsqlin(C, d, A, b, [], [], LB, UB, [], OptOptions);

Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Z_to_remove_ca = reshape(d, size(Z_to_remove_ca));

order = 6;
p = polyFit2D(T, X_P(:), Y_P(:), order, order);
y3 = polyVal2D(p,  X_P(:), Y_P(:), order, order);
y3(y3<0) = 0;

y3_residual = (Z_to_remove_ca(:) - C*y3);
y3_residual = y3_residual - nanmin(y3_residual);

LB_r = y3 - T;
LB_r(LB_r > 0) = 0;
[A, b] = DwellTime2D_CLLS_Assemble_A_b(m, n, 10);
T = lsqlin(C, y3_residual, A, b, [], [], LB_r, UB, [], OptOptions);
T = T + y3;

% Clear aperture results
z_removal_ca = C * T;
z_residual_ca = d - z_removal_ca;
z_removal_ca = reshape(z_removal_ca, size(Z_to_remove_ca));
z_residual_ca = reshape(z_residual_ca, size(Z_to_remove_ca));
Z_to_remove_ca = reshape(d, size(Z_to_remove_ca));

% Detilt
Z_to_remove_ca = RemoveSurface1(X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Z_to_remove_ca);
Z_to_remove_ca = Z_to_remove_ca - nanmin(Z_to_remove_ca(:));

z_removal_ca = RemoveSurface1(X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), z_removal_ca);
z_removal_ca = z_removal_ca - nanmin(z_removal_ca(:));

z_residual_ca = RemoveSurface1(X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e), z_residual_ca);


% Full results
z_removal = C_T * T;
z_residual = Z_to_remove(:) - z_removal;
z_removal = reshape(z_removal, size(Z_to_remove));
z_residual = reshape(z_residual, size(Z_to_remove));

% Reshape the T
T = reshape(T, size(X_P));

end