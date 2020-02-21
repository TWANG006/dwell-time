function Q = LegendreP2D_XYJC(x, y, j, c)
% Function
%   Q = LegendreP2D_XYJC(x, y, j, c)
% Purpose
%   Compute the 2D Orthoganal Legendre polynomials with order vctor J and its
%   corresponding coefficients vector C.
% Inputs: 
%   x, y: x, y coordinates of the Legendre polynomials normalized to [-1,1]
%      j: array of orders No. shown below
%      c: coefficients associated with each order j in J
% Outputs:
%      Q: The simulated surface error map
% Details:
%            1    n   ( n )       n-k        k
%   Pn(x) = --- Sigma (   )(x - 1)    (x + 1)
%           2^n  k=0  ( k )
%
%   Ln(y) = sqrt(2n+1)Pn(y), Lm(X) = sqrt(2m+1)Pn(y)
%
%   Qj(x,y) = Ln(x)Lm(y), Z = sum(CjQj(x,y))
%
%   No.j    n      m       Aberration Name              Polynomials  
%   Q1      0      0       piston                       1
%   Q2      1      0       x-tilt                       sqrt(3)x
%   Q3      0      1       y-tilt                       sqrt(3)y
%   Q4      2      0       x-power                      sqrt(5)/2(3x^2-1)
%   Q5      1      1       45 degree astigmatism        3xy
%   Q6      0      2       y-power                      sqrt(5)/2(3y^2-1)
%   Q7      3      0       x-primary coma               sqrt(7)/2(5x^3-3x)
%   Q8      2      1         
%   Q9      1      2
%   Q10     0      3       y-primary coma               sqrt(7)/2(5y^3-3y)
%   Q11     4      0       x-primary spherical          3/8(35x^4-30x^2+3)
%   Q12     3      1
%   Q13     2      2
%   Q14     1      3
%   Q15     0      4       y-primary spherical          3/8(35y^4-30y^2+3)
%
% Reference:
%   Mahajan, V. N. (2010). Orthonormal aberration polynomials for 
%   anamorphic optical imaging systems with rectangular pupils. 
%   Applied optics, 49(36), 6924-6929.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

%% 1. Deduce n, m from orders J
% 1    Q1(m=0,n=0)
% 2    Q2(m=1,n=0)  Q3(m=0,n=1) 
% 3    Q4(m=2,n=0)  Q5(m=1,n=1)  Q6(m=0,n=2)
% ...  Qj           Qj+1         Qj+2         ...
% ...
j = j(:);   % Column vector
c = c(:);   % Column vector

% Qj is in: row a & col b
a = floor((1 + sqrt(1 + 8*(j-1)))*0.5);
b = (j-1) - 0.5*(a.*(a-1)) + 1;

% Get n, m from a & b
m = a - b;
n = b - 1;

%% 2. Compute Legendre polynomial for each [n,m] pair
Q_partial = zeros(length(y),length(x), length(j));
for i=1:length(j)
    Q_partial(:,:,i) = LegendreP2D_XYMN(x, y, m(i), n(i));
    Q_partial(:,:,i) = Q_partial(:,:,i) * c(i);
end

Q = sum(Q_partial, 3);

end