function R = RotMat_handle()
% 3D rotation matrix for z-axis rotation.
    R = @(alpha) [cosd(alpha), -sind(alpha), 0 ;
                  sind(alpha), cosd(alpha), 0 ;
                  0, 0, 1;];
end
