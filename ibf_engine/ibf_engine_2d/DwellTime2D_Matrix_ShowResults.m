function DwellTime2D_Matrix_ShowResults(...
    X, Y,... the meshgird of the entire aperture [m]
    Z_to_remove,...height to remove in the full aperture [m]
    Z_removal,...the removal amount of the entire aperture [m]
    Z_residual,...the residual after IBF in the entire aperture [m]
    Z_to_remove_ca,...height to remove in the clear aperture [m]
    Z_removal_ca,...the removal amount in clear aperture [m]
    Z_residual_ca, ...the residual after IBF in clear aperture [m]
    X_P, Y_P,...the dwell positions [m]
    T, ...the dwell time map on the dwell positions [s]
    ca_range...the pixel size of the dwell grid [m]
    )
% Display the dwell time map defiend by X, Y, and Z with title

%% Copmute the sizes
% dx = BRF_params.d_pix;
% dy = BRF_params.d_pix;

rx = 0;
ry = 0;

%% dwell time
X_P_mm = X_P(ry+1:end-ry, rx+1:end-rx)*1e3;
Y_P_mm = Y_P(ry+1:end-ry, rx+1:end-rx)*1e3;
T = T(ry+1:end-ry, rx+1:end-rx);

%% Entire aperture to mm and nm units
X_mm = X*1e3;
Y_mm = Y*1e3;
Z_to_remove_nm = Z_to_remove*1e9;
Z_removal_nm = Z_removal*1e9;
Z_residual_nm = Z_residual*1e9;

%% Clear aperture to mm and nm units
X_ca_mm = X(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e)*1e3;
Y_ca_mm = Y(ca_range.y_s:ca_range.y_e, ca_range.x_s:ca_range.x_e)*1e3;
Z_to_remove_ca_nm = Z_to_remove_ca * 1e9;
Z_removal_ca_nm = Z_removal_ca*1e9;
Z_residual_ca_nm = Z_residual_ca*1e9;

%% Display the results
figure;

% Height to remove in the entire aperture
subplot(332);
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

% Height to remove in the clear aperture
subplot(333);
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

% Dwell time in on the dwell grid
subplot(334);
s = surf(X_P_mm, Y_P_mm, T);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Time [s]';
title({'Dwell Time Map',...
      ['Total Time = ' num2str(sum(T(:)) / 60.0) ' [min], ' ]},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

% Height removed in the entire aperture
subplot(335);
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

% Height removed in the clear aperture
subplot(336);
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

% Residual in the entire aperture
subplot(338);
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

% Residual in the clear aperture
subplot(339);
s = surf(X_ca_mm, Y_ca_mm, Z_residual_ca_nm);
s.EdgeColor = 'none';
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({'Height Residual in Clear Aperture',...
      ['PV = ' num2str(SmartPTV(Z_residual_ca_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_residual_ca_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;
view([0 90]);

end

