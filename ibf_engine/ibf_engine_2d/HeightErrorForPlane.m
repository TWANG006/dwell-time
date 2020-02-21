function Z_E = HeightErrorForPlane(Z)
   Z_mean = nanmean(Z(:));
   Z_E = Z - Z_mean;
   Z_E = Z_E - nanmean(Z_E(:));
end