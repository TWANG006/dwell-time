function Q = LegendreP2D_XYMN(x, y, m, n)
% Function:
%   Q = LegendreP2D_XYMN(x, y, m, n)
% Purpose:
%   Compute the 2D Orthoganal Legendre polynomials with order vector(m, n)
% Inputs: 
%   X, Y: x, y coordinates of the Legendre polynomials normalized to [-1,1]
%      m: order in x direction
%      n: order in y direction
% Outputs:
%      Q: Orthogonal Legendre surface
% Details:
%            1    n   ( n )       n-k        k
%   Pn(x) = --- Sigma (   )(x - 1)    (x + 1)
%           2^n  k=0  ( k )
%
%   Ln(y) = sqrt(2n+1)Pn(y), Lm(x) = sqrt(2m+1)Pn(x)
%
%   Qj(x,y) = Ll(X)Lm(y)
%
% Reference:
%   Mahajan, V. N. (2010). Orthonormal aberration polynomials for 
%   anamorphic optical imaging systems with rectangular pupils. 
%   Applied optics, 49(36), 6924-6929.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

Ln = legendreP(n, y);
Lm = legendreP(m, x);

Q = Ln(:)*Lm(:).';

end