function s = FWHM2Sigma(FWHM)
% Filname: 'FWHM2Sigma.m'. This function can be used for directly converting
% the full width at half maximum (FWHM) to the standard deviation of a peak.
%
% s = FWHM./(2*sqrt(2*log(2)))
% 
% One input arguments: 'FWHM'
% One output argument: 'sigma'
% 
% FWHM  :    Full width at half maximum
% s     :    Standard deviation
% 
% Input syntax: sigma = fwhm2sigma(FWHM)
% 
% Limitations
if nargin < 1 || nargin > 1
    error('Irregular amount of input arguments');
end
if FWHM <= 0,
    error('Full width at half maximum (FWHM) should be larger than 0');
end
% Calculation
s = FWHM./(2*sqrt(2*log(2)));
if s <= 0,
    error('Standard deviation (sigma) should be larger than 0');
end
end