function const = PrepConst
    const = struct;
    
    % Clayff potential
    % https://pubs.acs.org/doi/10.1021/jp0363287
    % Sequence (if applicable): Al, O, H, Si
    % Data from: octahedral aluminum, bridging oxygen, hydroxyl hydrogen
    % (charge only), tetrahedral silicon
    
    % Coulombic interaction
    % Vacuum permitivity
    const.epsil0 = 8.85419e-12; % Unit: F/m
    % Unit charge
    const.e = 1.602176634e-19; % Unit: Coulomb
    % Partial charge
    const.q = [1.5750, -1.0500, 0.4250, 2.100]; % Unit: e
    
    % vdW interaction
    % Emperical energy parameter, D_0 (Effective Hamaker constant)
    const.D = [1.3298e-6, 0.1554, 0, 1.8405e-6]; % Unit: kcal/mol
    % Emperial distance parameter, R_0
    const.R = [4.7943, 3.5532, 0, 3.7064]; % Unit: Angstrom
    
    const.D_convert = 6.9477e-21; % J = 1 kcal/mol
    
    % Polynomial coefficients for accelerating calculation
    const.coef12 = CalcPolyCoeffs(12);
    const.coef6 = CalcPolyCoeffs(6);
end