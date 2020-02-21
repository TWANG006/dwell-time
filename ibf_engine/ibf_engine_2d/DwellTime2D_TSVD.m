function [T ... [s] dwell map of the 2D IBF
    , z_removal... prediction of the removed height
    , z_residual ... prediction of the residual height
    , z_to_remove_ca...
    , z_removal_ca...predection of the removal amount using T
    , z_residual_ca...prediction of the residual after removal using T
    , C ... matrix C resembled from BRF
    ] = DwellTime2D_TSVD...
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

% This function implements the matrix-based Trancated SVD dwell time algorithm
%
%                 Nt
%   r(xk, yk) = Sigma C(xk - ui, yk-vi)T(ui, vi)
%                 i=1
%         k   ui'rd
%   T = sigma ------vi,  k <= nr, where si are the sigular values
%        i=1    si   

% 1. Dump X_P, Y_P dwell point positions into a 2D array as
%   |  u1    v1 |   P1
%   |  u2    v2 |   P2
%   | ...   ... |  ...
%   |  uNt   vNt|   PNt
P = [X_P(:), Y_P(:)];

% 2. Get the numbers of IBF machining points and sampling points of the surface error map R
Nt = size(P, 1);
Nr = numel(Z_to_remove);
% d_p = brfParams.d_pix;

% 3. Assemble the BRF matrix C, size(C) = Nr x Nt and vector d, Cx = d
[C, d, C_T] = DwellTime2D_Assemble_C_d(Nr, Nt, brfParams, Z_to_remove, X, Y, P, X_brf, Y_brf, Z_avg, ca_range, resampling_method);

% 4. SVD factorization of C
[U, S, V] = svd(C);
sigma = diag(S);            % Put the sigular values in a vector

% 5. TSVD algorithm
% delta = zeros(length(sigma), 1);
% rk = zeros(length(sigma), 1);

% Find the k_opt = min{k | rms(Z_residual_ca) <= rms_desired}
k_opt = 1;
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

for k = 1:length(sigma)
    tmp = U'*d;
    T = V*[tmp(1:k)./sigma(1:k); zeros(size(C,2)-k, 1)];
    
    Z_removal_ca = C * T;
    Z_residual_ca = Z_to_remove_ca(:) - Z_removal_ca;
    
    if(std(Z_residual_ca, 1) < rms_d)
        k_opt = k;
        break
    end
end
text = sprintf('Optimal Number of singular values = %d', k_opt);
disp(text);

% Piston adjustment for non-negativity
d_piston = d;

tmp = U'*d_piston;
T = V*[tmp(1:k_opt)./sigma(1:k_opt); zeros(size(C,2)-k_opt,1)];

while(min(T) < 0)
    d_piston = d_piston + 0.1e-9;
    tmp = U'*d_piston;
    T = V*[tmp(1:k_opt)./sigma(1:k_opt); zeros(size(C,2)-k_opt,1)];
end

text = sprintf('The minimum dwell time = %s', min(T));
disp(text);

% 6. Results
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