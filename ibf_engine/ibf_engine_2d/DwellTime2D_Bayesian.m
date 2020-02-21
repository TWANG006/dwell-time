function [B... BRF
    , X_B, Y_B...BRF Coordinates
    , X, Y...full aperture coordinates [m]
    , Z_removal, Z_residual... full aperture results [m]
    , T_P... dwell time on the dwell grid [s]
    , X_P, Y_P...dwell grid 
    , X_dw, Y_dw ... dwell grid coordinates [m]
    , Z_to_remove_dw, Z_removal_dw, Z_residual_dw...dwell grid results [m]
    , X_ca, Y_ca... clear aperture coordinates [m]
    , Z_to_remove_ca, Z_removal_ca, Z_residual_ca...[m]
    ] = DwellTime2D_Bayesian...
    ( Z_to_remove... height to remove [m]
    , BRF_mode...
    , BRF_params... BRF parameters
    , X_BRF, Y_BRF, Z_BRF ... averaged z & its coords
    , ca_range... Clear aperture range [pixel] 
    , pixel_m... pixel size [m/pixel]
    , options... optimization options
    )
% Purpose:
%   Implement the Bayesian-based dwell time algorithm proposed by Jiao et.
%   al in 2009
%
% Reference:
%   Jiao, C., Li, S., & Xie, X. (2009). Algorithm for ion beam figuring of 
%   low-gradient mirrors. Applied Optics, 48(21), 4090-4096.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

%% 0. Set the default options for the function
defaultOptions = struct(...
    'maxIters', 10,...
    'RMS_dif', 0.04e-9,...[m]
    'dwellTime_dif', 20,...[s]
    'lambda', 0.1e1,... % parameter for the total variation (TV) penalty
    'isDownSampling', false...
);
% 'alpha', 1e-1,... % step size for the additive iterations


%% 1. Deal with input arguments
% check invalid method
if nargin == 6
    options = defaultOptions; 
end
 
%% 2. Construct the BRF using the BRF parameters
% Release the BRF parameters
A = BRF_params.A;               % peak removal rate [m/s]
sigma_xy = BRF_params.sigma_xy; % standard deviation [m]
mu_xy = [0,0];                  % center is 0 [m]
d_p = BRF_params.d_pix;         % diameter [pixel]
r_p = floor(0.5 * d_p);         % radius [pixel]
l = options.lambda;             % parameter for the TV
% a = options.alpha;            % step size


% X, Y coordinates for B
[X_B, Y_B] = meshgrid(-r_p : r_p - 1, -r_p : r_p - 1);        
X_B = X_B * pixel_m;
Y_B = Y_B * pixel_m;

% Get B
if strcmpi(BRF_mode, 'avg')
    B = interp2(X_BRF, Y_BRF, Z_BRF, X_B, Y_B);
else
    B = BRFGaussian2D(X_B, Y_B, 1, [A, sigma_xy, mu_xy]); 
end 
%% 3. Define the dwell grid
% Get the size of the full aperture
[mM, nM] = size(Z_to_remove);   

% Get the dwell grid pixel range
dw_range.x_s = ca_range.x_s - r_p;   dw_range.x_e = ca_range.x_e + r_p;
dw_range.y_s = ca_range.y_s - r_p;   dw_range.y_e = ca_range.y_e + r_p;

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

%% 4. The real Bayesian algorithm
% The ca in dw range
% ca_in_dw_y_s = ca_range.y_s - dw_range.y_s;
% ca_in_dw_x_s = ca_range.x_s - dw_range.x_s;
% ca_in_dw_y_e = ca_in_dw_y_s + ca_range.y_e - ca_range.y_s;
% ca_in_dw_x_e = ca_in_dw_x_s + ca_range.x_e - ca_range.x_s;

% Get the Z_to_remove_dw in dwell grid and keep its non-negativity
Z_to_remove_dw = Z_to_remove(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
Z_to_remove_dw = Z_to_remove_dw - nanmin(Z_to_remove_dw(:));
Z_to_remove_ca = Z_to_remove(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
[x_ca,y_ca] = meshgrid(1:size(Z_to_remove_ca,2),1:size(Z_to_remove_ca, 1));
Z_to_remove_ca = RemoveSurface1(x_ca, y_ca, Z_to_remove_ca);
Z_to_remove_ca = Z_to_remove_ca - nanmin(Z_to_remove_ca(:));

% Bayesian algorithms with Total Variation regularizationm, assuming the
% noise model is Gaussian, i.e R(x) = E(x) * B(x) + NG(x), 
% a. Initialization
Tk0 =  Z_to_remove_dw;
Tk  = Tk0;
Bhat = B(end:-1:1, end:-1:1);
Bn = sum(B(:));

% b. Additive update
%   Tk+1 = Tk + a[(B`*E) - (B`*B)*Tk + ldiv(gradient_Tk / |gradient_Tk|)]
[GXk, GYk] = gradient(Tk);
norm_Gk = GXk.^2+GYk.^2;
norm_Gk = sqrt(sum(norm_Gk(:)));
GXk = GXk/norm_Gk;
GYk = GYk/norm_Gk;

% T_P = Tk./(l * divergence(GXk, GYk)).*(ConvFFT2(Tk, ConvFFT2(Bhat/Bn,B)) - ConvFFT2(Z_to_remove_dw, Bhat/Bn));
% T_P(T_P<0)=0;

% tmp = a*(ConvFFT2(Tk0, Bhat/Bn) - ConvFFT2(Tk, ConvFFT2(Bhat/Bn,B/Bn)) - l/a * (divergence(GXk, GYk)));
% T_P = Tk + tmp;
% T_P(T_P<0)=0;

% T_P = Tk + a * (-1 + l/a*divergence(GXk, GYk) + ConvFFT2(Tk0./ConvFFT2(Tk, B), Bhat/Bn));
% T_P(T_P<0)=0;

T_P = (Tk./(1 - l * divergence(GXk, GYk))) .* (ConvFFT2(Tk0./ConvFFT2(Tk, B), Bhat / Bn));


% figure;
% mesh(T_P);

% Dwell time results
T = zeros(size(Z_to_remove));;
T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_P;

% Calculate the height removal in the entire aperture
Z_removal = ConvFFT2(T,B); 
Z_residual = Z_to_remove - Z_removal;

% Obtain the dwell grid result
Z_removal_dw = Z_removal(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
Z_residual_dw = Z_residual(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

% Obtain the clear aperture results
Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);    
Z_residual_ca = Z_residual(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);

Tk = T_P;
Z_residual_ca_pre = Z_residual_ca;
num_iter = 1;

while(true)
    [GXk, GYk] = gradient(Tk);
    norm_Gk = GXk.^2+GYk.^2;
    norm_Gk = sqrt(sum(norm_Gk(:)));
    GXk = GXk/norm_Gk;
    GYk = GYk/norm_Gk;

%     T_P = Tk./(l * divergence(GXk, GYk)).*(ConvFFT2(Tk, ConvFFT2(Bhat/Bn,B)) - ConvFFT2(Z_to_remove_dw, Bhat/Bn));
%     T_P(T_P<0)=0;
    
%     tmp = a*(ConvFFT2(Tk0, Bhat/Bn) - ConvFFT2(Tk, ConvFFT2(Bhat/Bn,B/Bn)) - l/a * (divergence(GXk, GYk)));
%     T_P = Tk + tmp;
%     T_P(T_P<0)=0;
    
%     T_P = Tk + a * (-1 + l/a*divergence(GXk, GYk) + ConvFFT2(Tk0./ConvFFT2(Tk, B), Bhat/Bn));
%     T_P(T_P<0)=0;

    T_P = (Tk./(1 - l * divergence(GXk, GYk))) .* (ConvFFT2(Tk0./ConvFFT2(Tk, B), Bhat / Bn));
    
%     mesh(T_P);
%     drawnow;
%     display(sum(T_P(:))/60);

    % Dwell time results
    T = zeros(size(Z_to_remove));
    T(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e) = T_P;

    % Calculate the height removal in the entire aperture
    Z_removal = ConvFFT2(T,B); 
    Z_residual = Z_to_remove - Z_removal;

    % Obtain the dwell grid result
    Z_removal_dw = Z_removal(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);
    Z_residual_dw = Z_residual(dw_range.y_s:dw_range.y_e, dw_range.x_s:dw_range.x_e);

    % Obtain the clear aperture results
    Z_removal_ca = Z_removal(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);    
    Z_residual_ca = Z_residual(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e);
    
    if(nanstd(Z_residual_ca(:) - Z_residual_ca_pre(:), 1) < options.RMS_dif)
        display('RMS difference limit reached.');
        display(['Number of iterations = ', num2str(num_iter)]);
        break;
    elseif(nanstd(Z_residual_ca(:),1) > nanstd(Z_residual_ca_pre(:),1))
        display('RMS of the residual is getting worse, return the previous result.');
        display(['Number of iterations = ', num2str(num_iter)]);
        T_P = Tk;
        break;
    elseif(num_iter > options.maxIters)
        display('Maximum number of iterations reached.');
        break;
    elseif (abs(sum(T_P(:) - Tk(:))) / sum(Tk(:)) < options.dwellTime_dif)
        display('Threshold hitted');
        display(['Number of iterations = ', num2str(num_iter)]);
        break;
    end    
    
    Tk = T_P;
    Z_residual_ca_pre = Z_residual_ca;
    num_iter = num_iter +1;
end

Z_removal_ca = RemoveSurface1(x_ca, y_ca, Z_removal_ca);
Z_removal_ca = Z_removal_ca - nanmin(Z_removal_ca(:));
Z_residual_ca = RemoveSurface1(x_ca, y_ca, Z_residual_ca);


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
    [C, d, C_T] = DwellTime2D_FFT_Assemble_C_d(Nr, Nt, BRF_params, Z_to_remove, X, Y, P, ca_range);

    % Down sample T_dw
    T_P = imresize(T_P, 1/interval_P_m, 'bicubic') * interval_P_m.^2;
    T_P_v = T_P(:);
    T_P_v(T_P_v<0) = 0;

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
end

end