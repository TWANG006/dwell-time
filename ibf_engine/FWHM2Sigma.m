function s = FWHM2Sigma(FWHM)
% Function: 
%   s = FWHM2Sigma(FWHM)
%
% Purpose:
%   This function can be used for directly converting the full width at 
%   half maximum (FWHM) to the standard deviation of a peak.
%
%   s = FWHM./(2*sqrt(2*log(2)))
% 
% Input:
%   FWHM: Full width at half maximum
% 
% Output: 
%   s: Standard deviation
% 
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

s = FWHM./(2*sqrt(2*log(2)));

end