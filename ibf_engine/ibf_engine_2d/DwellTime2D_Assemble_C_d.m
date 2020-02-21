function [C, z, C_T] = DwellTime2D_Assemble_C_d(...
    Nr, ... number of elements in Z_to_remove
    Nt, ... number of dwell positions
    BRF_params,... BRF parameters
    Z_to_remove,... desired height to remove
    X, Y, ... X, Y coordinates of Z_to_remove
    P, ...  Dwell time positions
    X_brf, Y_brf, ... BRF X, Y coordinates
    Z_avg,... Averaged BRF 
    ca_range,...range of the clear aperture
    resampling_method... use 'model' or 'avg'
)
% Purpose:
%   Assemble the matrix C and vector d
%
% Matrix C is assembled as below:
%
%   | r1 |     | c11    c12    ...    c1Nt | | t1 |
%   | r2 |     | c21    c22    ...    c2Nt | | t2 |
%   | .  |     |  .      .      .       .  | | .  |
%   | .  |     |  .      .      .       .  | | .  |
%   | .  |     |  .      .      .       .  | | .  |
%   | rNr|     | c21    c22    ...    c2Nt | | tNt|
%
% Reference:
%   Zhou, L., Dai, Y. F., Xie, X. H., Jiao, C. J., & Li, S. Y. (2007). 
%   Model and method to determine dwell time in ion beam figuring.
%   Nanotechnol. Precis. Eng., 5(8–9), 107-112.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

% 1. Release the BRF parameters
A = BRF_params.A;                   % peak removal rate [m/s]
sigma_xy = BRF_params.sigma_xy;     % standard deviation [m]
mu_xy = [0,0];                      % center is 0 [m]

% 2. Get the clear aperture size
ca_m = ca_range.y_e - ca_range.y_s + 1; 
ca_n = ca_range.x_e - ca_range.x_s + 1;

row_C= ca_m * ca_n;

% 3. Assemble the matirx C
C = zeros(row_C, Nt);
C_T = zeros(Nr, Nt);
for i = 1:Nt
    Yk = Y - P(i, 2);   % yk - vi
    Xk = X - P(i, 1);   % xk - ui
    
    if strcmpi(resampling_method, 'avg')
        z_brf = interp2(X_brf, Y_brf, Z_avg, Xk, Yk);
        z_brf(~isfinite(z_brf))=0;
    elseif strcmpi(resampling_method, 'model')
        z_brf = BRFGaussian2D(Xk, Yk, 1,[A, sigma_xy, mu_xy]);
    end
    
    C_T(:, i) = z_brf(:);
    z_brf = z_brf(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
    C(:, i) = z_brf(:);
end    

% 4. Assemble the vector d
z = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
z = z - nanmin(z(:));
z = z(:);

end
