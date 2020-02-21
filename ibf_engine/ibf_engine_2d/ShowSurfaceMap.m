function ShowSurfaceMap(X, Y, Z, title_str)
% Function:
%   ShowSurfaceMap(X, Y, Z, title_str)
%
% Purpose:
%   Display the surface map defiend by X, Y, and Z with title
%
% Inputs:
%        X, Y: meshgrid points of the surface map [m]
%           Z: the surface map [m]
%   title_str: the title string
%
% Info:
%   Contact: tianyiwang666@gmail.com (Dr WANG Tianyi)
%   Copyright reserved.

% change to mm and nm units
X_mm = X*1e3;
Y_mm = Y*1e3;
Z_nm = Z*1e9;

% display the map
figure;
pcolor(X_mm, Y_mm, Z_nm), shading interp;
axis image; 
c = colorbar;
c.Label.String = 'Height [nm]';
title({title_str,...
      ['PV = ' num2str(SmartPTV(Z_nm(:))) ' nm, ' 'RMS = ' num2str(nanstd(Z_nm(:),1)) ' nm']},...
    'FontWeight', 'normal',...
    'FontSize', 12);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap jet;

end