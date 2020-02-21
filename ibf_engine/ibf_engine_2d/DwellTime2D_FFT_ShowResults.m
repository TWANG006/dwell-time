function DwellTime2D_FFT_ShowResults(B,... BRF
    X_B, Y_B,...BRF Coordinates
    X, Y,... the meshgird of the entire aperture [m]
    Z_to_remove,... the amount to remove in the entire aperture [m]
    Z_removal,... the removal amount of the entire aperture [m]
    Z_residual,... the residual after IBF in the entire aperture [m]
    T_P, ...the dwell time map on the dwell positions [s]
    X_P, Y_P,...
    X_dw, Y_dw,...the dwell grid [m]
    Z_to_remove_dw,... the amount to remove in the entire aperture [m]
    Z_removal_dw,... the removal amount of the entire aperture [m]
    Z_residual_dw,... the residual after IBF in the entire aperture [m]
    X_ca, Y_ca,... the dwell positions [m]
    Z_to_remove_ca,...the amount to remove in clear aperture [m]
    Z_removal_ca,...the removal amount in clear aperture [m]
    Z_residual_ca ...the residual after IBF in clear aperture [m]
    )
% Display the final results of the FFT-based dwell time algorithm

%% Units change
% BRF
X_B_mm = X_B * 1e3;
Y_B_mm = Y_B * 1e3;
B_nm = B * 1e9;

% Entire aperture
X_mm = X * 1e3;
Y_mm = Y * 1e3;
Z_to_remove_nm = Z_to_remove * 1e9;
Z_removal_nm = Z_removal * 1e9;
Z_residual_nm = Z_residual * 1e9;

% Clear aperture
X_ca_mm = X_ca * 1e3;
Y_ca_mm = Y_ca * 1e3;
Z_to_remove_ca_nm = Z_to_remove_ca * 1e9;
Z_removal_ca_nm = Z_removal_ca * 1e9;
Z_residual_ca_nm = Z_residual_ca * 1e9;

% Dwell time
X_P_mm = X_P * 1e3;
Y_P_mm = Y_P * 1e3;

% Dwell grid
X_dw_mm = X_dw * 1e3;
Y_dw_mm = Y_dw * 1e3;
Z_to_remove_dw_nm = Z_to_remove_dw * 1e9;
Z_removal_dw_nm = Z_removal_dw * 1e9;
Z_residual_dw_nm = Z_residual_dw * 1e9;

%% Display the maps
fsfig('Dwell time calculation results');

%% BRF
subplot(3,4,1);
s = surf(X_B_mm, Y_B_mm, B_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title('BRF',...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

%% Dwell Time 
% 2. Dwell time Tdw for the dwell grid
subplot(3,4,5);
s = surf(X_P_mm, Y_P_mm, T_P);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Time [s]';
title({'Dwell Time on Dwell Grid',...
      ['Total Time = ' num2str(sum(T_P(:)) / 60.0) ' [min], ' ]},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

%% Entire aperture
% 4. Height to remove in the entire aperture
subplot(3,4,2);
s = surf(X_mm, Y_mm, Z_to_remove_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height to remove',...
      ['PV = ' num2str(SmartPTV(Z_to_remove_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_to_remove_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 5. Height removed in the entire aperture
subplot(3,4,6);
s = surf(X_mm, Y_mm, Z_removal_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Removed',...
      ['PV = ' num2str(SmartPTV(Z_removal_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_removal_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 6. Residual in the entire aperture
subplot(3,4,10);
s = surf(X_mm, Y_mm, Z_residual_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Residual',...
      ['PV = ' num2str(SmartPTV(Z_residual_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_residual_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

%% Dwell Grid
% 7. Height to remove in the dwell grid
subplot(3,4,3);
s = surf(X_dw_mm, Y_dw_mm, Z_to_remove_dw_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height to remove in dwell grid',...
      ['PV = ' num2str(SmartPTV(Z_to_remove_dw_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_to_remove_dw_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 8. Height removed in the dwell grid
subplot(3,4,7);
s = surf(X_dw_mm, Y_dw_mm, Z_removal_dw_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Removed in dwell grid',...
      ['PV = ' num2str(SmartPTV(Z_removal_dw_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_removal_dw_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 9. Residual in the dwell grid
subplot(3,4,11);
s = pcolor(X_dw_mm, Y_dw_mm, Z_residual_dw_nm);
[mca, nca] = size(X_ca_mm);
pixel_mm = X_ca_mm(1, 2) - X_ca_mm(1, 1); 
rectangle('Position', [X_ca_mm(1,1), Y_ca_mm(1,1), nca * pixel_mm, mca * pixel_mm]);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Residual in dwell grid',...
      ['PV = ' num2str(SmartPTV(Z_residual_dw_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_residual_dw_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

%% Clear Aperture
% 10. Height to remove in the clear aperture
subplot(3,4,4);
s = surf(X_ca_mm, Y_ca_mm, Z_to_remove_ca_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height to remove in clear aperture',...
      ['PV = ' num2str(SmartPTV(Z_to_remove_ca_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_to_remove_ca_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 11. Height removed in the clear aperture
subplot(3,4,8);
s = surf(X_ca_mm, Y_ca_mm, Z_removal_ca_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Removed in Clear Aperture',...
      ['PV = ' num2str(SmartPTV(Z_removal_ca_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_removal_ca_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% 12. Residual in the clear aperture
subplot(3,4,12);
s = surf(X_ca_mm, Y_ca_mm, Z_residual_ca_nm);
s.EdgeColor = 'none';
axis image; 
rms_Z = nanstd(Z_residual_ca_nm(:),1);
caxis([-1 1]*3*rms_Z);
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Residual in Clear Aperture',...
      ['PV = ' num2str(SmartPTV(Z_residual_ca_nm(:))) ' nm, ' 'RMS = ' num2str(rms_Z) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

end