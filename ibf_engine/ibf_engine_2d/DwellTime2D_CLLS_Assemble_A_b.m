function [A, b] = DwellTime2D_CLLS_Assemble_A_b(m, n, max_dt_diff)
% Purpose:
%   Assemble the constraints matrix A, the size of the dwell points 
%   map is = m x n, size(A) = [(m-1)n+m(n-1)] x Nt. Suppose m = 4, n = 3, 
%   then Nt = 12, the dwell point positions are like: 
%       | P1  P5   P9 |
%       | P2  P6  P10 |
%       | P3  P7  P11 |
%       | P4  P8  P12 |
%   Then, the first half of A is assembled by two kinds of BLOCKS 
%   seperated by "=" like:
% =======================================            ===========
% | -1  1                               |            | P2 - P1  |
% |    -1  1                            |            | P3 - P2  |  
% |       -1  1                         | | P1  |    | P4 - P3  |
% --------------------------------------- | P2  |    -----------
% |             -1  1                   | | P3  |    | P6 - P5  | 
% |                -1  1                | | P4  |    | P7 - P6  |
% |                   -1  1             | | P5  |    | P8 - P7  |
% --------------------------------------- | P6  |    -----------
% |                         -1  1       | | P7  | =  | P10 - P9 |
% |                            -1  1    | | P8  |    | P11 - P10|
% |                               -1  1 | | P9  |    | P12 - P11|
% ======================================= | P10 |    ============
% |-1            1                      | | P11 |    | P5 - P1  |
% |   -1            1                   | | P12 |    | P6 - P2  |
% |      -1            1                |            | P7 - P3  |
% |         -1            1             |            | P8 - P4  |
% ---------------------------------------            ------------
% |             -1           1          |            | P9 - P5  |
% |                -1           1       |            | P10 - P6 |
% |                   -1           1    |            | P11 - P7 |
% |                      -1           1 |            | P12 - P8 |
% =======================================            ============
% Second half can just be done as
%       [  A  ]
% A =   [     ]
%       [ -A  ]
% Reference:
%   Wang, T., Huang, L., Vescovi, M., Kuhne, D., Tayabaly, K., 
%   Bouet, N., & Idir, M. (2019). Study on an effective one-dimensional 
%   ion-beam figuring method. Optics express, 27(11), 15368-15381.
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.
size_p = m * n;
exp_p = [zeros(m,n);NaN(1, n)]; % Expand matrix p to add NaN at bottom

% Put the -1 & 1 along the diagonal
e = ones(size_p, 1);
Dx = spdiags([-e, e], [0, m], m*(n-1), size_p);
Dy = spdiags([-e, e], [0, 1], size_p, size_p);

Gy = exp_p(2:end, :) - exp_p(1:end-1, :);

% Construct A half-by-half
A = [Dy(isfinite(Gy),:);Dx];
A = [A;-A];
A = full(A);
% Assemble b
b = ones(size(A,1), 1)*max_dt_diff;

end