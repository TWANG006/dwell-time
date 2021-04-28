function FWHM = Sigma2FWHM(sigma)
% Function: 
%   FWHM = Sigma2FWHM(sigma)
%
% Purpose:
%   This function can be used for directly converting the standard 
%   deviation of a peak to the full width at half maximum (FWHM).
% 
% Input:
%   s: Standard deviation
% 
% Output: 
%   FWHM: Full width at half maximum
% 
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

FWHM = 2*sqrt(2*log(2))*sigma;

end