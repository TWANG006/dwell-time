function smartPTV = SmartPTV(x)
% Function:
%   smartPTV = SmartPTV(x)
%
% Input:
%   x: the array that need to calculate PV
%
% Output:
%   smartPTV: the scalar PV value of x
%
% Purpose:
%   Compute PV value without considering any NaNs
%
% Info:
%   Contact: huanglei0114@gmail.com (Dr. HUANG Lei), tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

smartPTV = range(x(isfinite(x)));