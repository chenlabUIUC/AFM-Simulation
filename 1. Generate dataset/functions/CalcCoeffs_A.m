function coef = CalcCoeffs_A(const, bias, pos, atom_pos, cutoff_range, convert)
    coef = zeros(1,21);
    
    select = atom_pos(:,2:3) - pos(1:2);
    select = sqrt(sum(select.^2,2)) < cutoff_range;
    pos2 = atom_pos(select,:); % list for neighboring atoms
    pos2(:,2:4) = pos2(:,2:4) - pos; % xy relative position in px
%     pos2(:,2:4) = pos2(:,2:4) * convert; % now in Å
    
    q = const.q + bias(:,1)';
    D = const.D + bias(:,2)';
    D = D * const.D_convert; % J
    R = const.R + bias(:,3)';    
    R = R * 1e-10; % m
    
    % Can use matrix calculation.
        r = sqrt(sum(pos2(:,2:4).^2,2)) * 1e-10; % m
        dz = pos2(:,4) * 1e-10; % m
        qj = q(pos2(:,1))';
        Dj = D(pos2(:,1))';
        Rj = R(pos2(:,1))';
        
        qi = 2.100;
        Di = 1.8405e-6 * const.D_convert;
        Ri = 3.7064;
    % Coulomb
        temp = - dz ./ (r.^3) .* qj * const.e^2 / 4 / pi / const.epsil0;
        coef(1) = sum(temp);
%         coef(1) = sum(temp) * qi;
    % vdW 12
        for i = 1:13
            temp = const.coef12(i) .* Rj.^(i-1) .* 1e-10^(13-i)...
                   .* sqrt(Dj) .* dz ./ (r.^14) * -12;
            coef(1+i) = sum(temp);
%             coef(1+i) = sum(temp) * sqrt(Di) * Ri^(13-i);
        end
    % vdW 6
        for i = 1:7
            temp = const.coef6(i) .* Rj.^(i-1) .* 1e-10^(7-i)...
                   .* sqrt(Dj) .* dz ./ (r.^8) * -6;
            coef(14+i) = sum(temp);
%             coef(14+i) = sum(temp) * sqrt(Di) * Ri^(7-i);
        end
end