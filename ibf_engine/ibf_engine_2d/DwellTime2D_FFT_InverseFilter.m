function T = DwellTime2D_FFT_InverseFilter(R, B, gamma, useDCT)
% 2D Inverse filtering for deconvolution
% Inputs:
%   R: the filtered signal
%   B: the kernel
%   n: the filtering threshold
% Output:
%   T: the recovered original signal

%% Size of the Signal and the Kernel
[mR, nR] = size(R);

%% Padding and FFT
if useDCT == true
    % Mirror the matrix
    R_top = [R, fliplr(R)];
    R_T  = [R_top; flipud(R_top)];
    
    % Compute the FFT padding size
    [mRt, nRt] = size(R_T);

    FR = fft2(R_T);
    FB = fft2(B, mRt, nRt);
else
    % Perform FFT
    FR = fft2(R);
    FB = fft2(B, mR, nR);
end

%% Thresholding
sFB = FB.*(abs(FB)>0) + 1/gamma*(abs(FB)==0);
iFB = 1./sFB;
iFB = iFB.*(abs(sFB)*gamma>1)+gamma*(abs(sFB)*gamma<=1);

%% Inverse filtering
T = real(ifft2(iFB.*FR));

if useDCT == true
    T = T(1:mR, 1:nR);
% else
%     T = ConvFFT2_Crop(T, mR, nR, mB, nB);    
end

%% Non-negative processing
% min_t = nanmin(T(:));
% if min_t <0
%     T = T + abs(min_t);
% end
T(T<0) = 0;

end