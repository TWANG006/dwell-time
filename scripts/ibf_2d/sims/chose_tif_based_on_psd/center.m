function [cen] = center(len)
    %%
    % calculates the center of an input number for array/matrix use.
    % Useful because fft's use this position as the center position.
    
    cen = ceil( (len + 1)/2 );