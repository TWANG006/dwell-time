function [T ... [s] dwell map of the 2D IBF
    , z_removal... prediction of the removed height
    , z_residual ... prediction of the residual height
    , z_to_remove_ca...
    , z_removal_ca...predection of the removal amount using T
    , z_residual_ca...prediction of the residual after removal using T
    , C ... matrix C resembled from BRF
    ] = DwellTime2D_LSQR...
    ( Z_to_remove... desired removal map
    , X... X grid points of the z_to_remove
    , Y... Y grid points of the z_to_remove
    , brfParams... BRF parameters
    , X_brf, Y_brf, Z_avg...
    , X_P... x coordinates of the dwell points, size = m x n
    , Y_P... y coordinates of the dwell points, size = m x n
    , ca_range...
    , rms_d... desired rms value of the final residual
    , resampling_method...
    )
% Purpose:
%   This function implements the matrix-based LSQR dwell time algorithm
%
%                 Nt
%   r(xk, yk) = Sigma C(xk - ui, yk-vi)T(ui, vi)
%                 i=1
%         k   ui'rd
%   T = sigma ------vi,  k <= nr, where si are the sigular values
%        i=1    si   
%
% Reference:
%   [1] Carnal, C. L., Egert, C. M., & Hylton, K. W. (1992, December). 
%       Advanced matrix-based algorithm for ion-beam milling of 
%       optical components. In Current Developments in Optical Design and 
%       Optical Engineering II (Vol. 1752, pp. 54-62). 
%       International Society for Optics and Photonics.
%   [2] Wu, J. F., Lu, Z. W., Zhang, H. X., & Wang, T. S. (2009). 
%       Dwell time algorithm in ion beam figuring. 
%       Applied optics, 48(20), 3930-3937.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

% 1. Dump X_P, Y_P dwell point positions into a 2D array as
%   |  u1    v1 |   P1
%   |  u2    v2 |   P2
%   | ...   ... |  ...
%   |  uNt   vNt|   PNt
P = [X_P(:), Y_P(:)];

% 2. Get the numbers of IBF machining points and sampling points of the surface error map R
Nt = size(P, 1);
Nr = numel(Z_to_remove);

% 3. Assemble the BRF matrix C, size(C) = Nr x Nt and vector d, Cx = d
[C, d, C_T] = DwellTime2D_Assemble_C_d(Nr, Nt, brfParams, Z_to_remove, X, Y, P, X_brf, Y_brf, Z_avg, ca_range, resampling_method);


% 4. Find the optimal damp factor fo LSQR algorithm
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
[m,n] = size(C);

% 5. Piston adjustment & damp optimization
gamma = 0;
damp = 1e-9;
d_piston = d + gamma;
[T,~,~,~,~,~,~,~,~,~] = lsqrSOL(m, n, C, d_piston, damp,0, 0, 0, 100,0);

while(min(T) < 0)
    % Piston adjuestment
    gamma = gamma + 0.1e-9;
    d_piston = d + gamma;
    
    % Optimize damp factor
    [T,~,~,~,~,~,~,~,~,~] = lsqrSOL(m, n, C, d_piston, damp,0, 0, 0, 100,0);
    z_removal_ca = C * T;
    z_residual_ca = Z_to_remove_ca(:) - z_removal_ca;

    while(std(z_residual_ca, 1) > rms_d)
    %     disp(std(z_residual_ca, 1));
        damp = damp - 0.01e-9;
        [T,~,~,~,~,~,~,~,~,~] = lsqrSOL(m, n, C, d_piston, damp,0, 0, 0, 100,0);
        z_removal_ca = C * T;
        z_residual_ca = Z_to_remove_ca(:) - z_removal_ca;
    end
    text = sprintf('Optimized damp factor = %e', damp);
    disp(text);
%     disp(min(T));    
%     [T,~,~,~,~,~,~,~,~,~] = lsqrSOL(m, n, C, d_piston, damp,0, 0, 0, 100,0);
end
% 
% disp(sum(T)/60);
text = sprintf('Piston adjustment done. The piston added is is %e [nm]', gamma*1e9);
disp(text);

% 6. Add path and surface error weights
% TODO


% 7. Results
% Clear aperture results
z_removal_ca = C * T;
z_residual_ca = d_piston - z_removal_ca;
z_removal_ca = reshape(z_removal_ca, size(Z_to_remove_ca));
z_residual_ca = reshape(z_residual_ca, size(Z_to_remove_ca));
z_to_remove_ca = reshape(d_piston, size(Z_to_remove_ca));

% Detilt
X_ca = X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
Y_ca = Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
z_to_remove_ca = RemoveSurface1(X_ca, Y_ca, z_to_remove_ca);
z_to_remove_ca = z_to_remove_ca - nanmin(z_to_remove_ca(:));

z_removal_ca = RemoveSurface1(X_ca, Y_ca, z_removal_ca);
z_removal_ca = z_removal_ca - nanmin(z_removal_ca(:));

z_residual_ca = RemoveSurface1(X_ca, Y_ca, z_residual_ca);

% Full results
z_removal = C_T * T;
z_residual = Z_to_remove(:) - z_removal;
z_removal = reshape(z_removal, size(Z_to_remove));
z_residual = reshape(z_residual, size(Z_to_remove));

% Reshape the T
T = reshape(T, size(X_P));

end