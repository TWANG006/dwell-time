function uk = SurfaceExtension_GerchbergPapoulis(...
    u0, ...initial extended signal
    G, ...spatial gate function
    Gy, ...spatial gate function in y direction only
    Gox, ...frequency-domain gate function in horizontal direction
    Goy, ...frequency-domain gate function in vertical direction
    rms_thrd,...the rms threshold between each two consecutive iterations
    max_iter...max number of iterations
    )
% Function to perform the improved 2D Gerchberg-Papoulis bandlimited
% surface extrapolation algorithm. 
%
% Refrence:
%   Marks, R. J. (1981). Gerchberg’s extrapolation algorithm in two 
%   dimensions. Applied optics, 20(10), 1815-1820. 

%% Initialization
if nargin == 5
    rms_thrd = 1e-9;
    max_iter = 500;
end
if nargin == 6
    max_iter = 500;
end

u_pre = u0.*G; v_pre = u_pre; w_pre = u_pre;

%% Iterative update
i=1;
while i <= max_iter
    wk = (1-G).*ifft(ifftshift(Gox.*fftshift(fft(w_pre, [], 2), 2), 2), [], 2);
    vk = (1-Gy).*ifft(ifftshift(Goy.*fftshift(fft(v_pre, [], 1), 1), 1), [], 1) + wk;
    uk = u_pre + vk;
    
    % Early stop when the rms difference is satisfied 
    if(nanstd(uk(:) - u_pre(:), 1) <= rms_thrd)
        break;
    end
   
    u_pre = uk;
    v_pre = vk;
    w_pre = wk;
    i=i+1;
end

end