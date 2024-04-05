%% ForceCalculation_calcMatrix.m
% Calculate the force matrix for fast calculation of tip-surface
% interaction afterwards. This step is very time-comsuming, so we generate
% 128*128 px images instead of 256*256 px, with buffer size.
% 
% INPUT:
%   'in_path': defect structures generated from 'AtomicSimulation_defect.m'
%   'convert': pixel size, Å/px
%   'plot_size_px': targeted image size to generate.
%   'plot_buffer_px': image buffer size.
%   'Bias': a handle to easily tuning the interaction parameters for atoms.
%           We set zeros here so just use the original Clayff potential.
%   'cutoff_range_A': cutoff distance for force calculation.
%   'Fz_init': init value for zero force point finding.
%   'Tip_radius': estimated tip size for plotting, not entering the
%                 calculation.
%   'Tip_height': estimated tip size for plotting, not entering the
%                 calculation.
%   'd_list': the discrete z height we calculated to get the force curve.
% 
% OUTPUT:
%   'F_mat': 4D matrix with size (n_xy, n_xy, n_d, 21), 
%       n_xy is the number of steps calculated in the x,y directions,
%       related to image size and buffer size.
%       n_d is the length of d_list.
%       21 is the number of parameters required to calculate the force at
%       each point following Clayff potential.
%   'd_list': same as input, saved for cleaner interpretation of F_mat.
% 
% NOTE: Atom type 1--Al, 2--O, 3--H, 4--Si
%
% HISTORY:  Written by Chang (Chris) Qian Last modified by Chang (Chris)
% Qian on 04/04/2024

clear; close all; clc;
addpath('functions/')

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool(6);
    p = gcp('nocreate');
end
pctRunOnAll warning('off','all')

in_path = 'output/Defected structure/';
out_path = 'output/Force matrix/';
flist = dir(in_path);
flist = flist(3:end);

for i = length(flist):-1:1

load([in_path, flist(i).name])

%%%% workflow control %%%%
check_configuration = 1;

save_mat_name = [out_path, flist(i).name];
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_size_px = 128; % px
plot_buffer_px = 10; % px
Const = PrepConst;

% for manually change interaction parameter
% Columns: Al, O, H, Si
% Rows: partial charge, energy parameter, distance parameter
Bias = [0,  0,  0,  0;
        0,  0,  0,  0;
        0,  0,  0,  0]';
    
cutoff_range_A = 10; % A
Fz_init = 3; % Å
Tip_radius = 5; % A
Tip_height = 6; % A

d_list = [-3:0.1:17]; % Å
convert = 200/256; % Å/px

%% Larger plot size for further calculation
plot_size_px = (plot_size_px*convert + (Tip_radius + cutoff_range_A)*2) / convert;
plot_size_px = round(plot_size_px) + plot_buffer_px * 2;

%% Read in, parameter initiation

% atom position (A)
atom_pos = output.atom_pos;
atom_size = output.atom_size;

% Probe x,y list
xy_1d = [1:plot_size_px] - 0.5;
xy_1d = xy_1d * convert; % A
[x,y] = meshgrid(xy_1d,xy_1d);
xy_2d = [reshape(x,[],1),reshape(y,[],1)]; % A

% interaction calculation
F_mat = zeros(length(d_list), plot_size_px * plot_size_px);
I0 = zeros(plot_size_px * plot_size_px,1);

% Generate tip
% Parameters: radius, height, lattice_points_per_side, rot111
Tip = generate_Si(Tip_radius,Tip_height,ceil(Tip_radius),1);

n_xy = length(xy_1d);
n_d = length(d_list);

%% Quick check on configuration:
if check_configuration
    figure(1); hold on;
    axis equal;
    xlim([min(atom_pos(:,2)),max(atom_pos(:,2))])
    ylim([min(atom_pos(:,3)),max(atom_pos(:,3))])

    % Atoms
    scatter3(atom_pos(:,2),atom_pos(:,3),atom_pos(:,4),...
             atom_size(atom_pos(:,1)).^2*30, atom_pos(:,4),'filled')

    % Generating area
    % scatter(xy_2d(:,1),xy_2d(:,2),'x')

    % Cutoff range
    % [X,Y,Z] = sphere;
    % surf(X.* cutoff_range_A,...
    %      Y.* cutoff_range_A,...
    %      Z.* cutoff_range_A)

    % Tip geometry
    Tip(:,:) = Tip(:,:) + [50,50,10];
    scatter3(Tip(:, 1), Tip(:, 2), Tip(:, 3), 'bo', 'filled');
        % Plot bonds between neighboring atoms on the tip
        lattice_constant = 5.43;
        for i = 1:size(Tip,1)
            for j = i+1:size(Tip,1)
                distance = norm(Tip(i, :) - Tip(j, :));
                if distance <= 1.001*sqrt(3)/4 * lattice_constant  % You can adjust the bond length threshold as needed
                    plot3([Tip(i, 1), Tip(j, 1)], ...
                          [Tip(i, 2), Tip(j, 2)], ...
                          [Tip(i, 3), Tip(j, 3)], 'k-');
                end
            end
        end
end

%% Generate matrix for further force calculation
% interaction calculation
F_temp = cell((n_xy * n_xy * n_d),1);

t3 = n_xy * n_xy * n_d / p.NumWorkers;
tic
parfor idx = 1:(n_xy * n_xy * n_d)
    if ~mod(idx,1e5)
        tic
    end
    
    [x_i, y_i, d_i] = ind2sub([n_xy, n_xy, n_d], idx);
    x = xy_1d(x_i);
    y = xy_1d(y_i);
    d = d_list(d_i);
    pos = [x, y, d];
    
    coefs = CalcCoeffs_A(Const, Bias, pos, atom_pos, cutoff_range_A, convert);
    
    F_temp{idx} = coefs;
    
    if ~mod(idx,1e5)
        t = toc;
        t2 = t * t3;
        disp([num2str(t,'%.3e'), 's  ', 'Estimated: ',...
              num2str(t2, '%.3e'), 's'])
    end
end

totalSeconds = toc;
hours = floor(totalSeconds / 3600);
minutes = floor((totalSeconds - hours * 3600) / 60);
seconds = mod(totalSeconds, 60);
formattedTime = datestr(datenum(0,0,0,hours,minutes,seconds), 'HH:MM:SS');
disp(['Total time: ', formattedTime])


F_mat = zeros(n_xy, n_xy, n_d, 21);
for x_i = 1:n_xy
    for y_i = 1:n_xy
        for d_i = 1:n_d
            idx = sub2ind([n_xy, n_xy, n_d], x_i, y_i, d_i);
            F_mat(x_i, y_i, d_i, :) = F_temp{idx};
        end
    end
end

%% Save result
save(save_mat_name,'F_mat','d_list')

end

%% Functions
function E = Energy_Coulombic(const, bias, pos_i, pos_j, type_i, type_j)
%     E = 0;
    q = const.q + bias(:,1)';
    
    r = pos_i-pos_j;
    r = sqrt(sum(r.^2)); % A
    r = r * 1e-10; % m
    
    E = const.e^2 / 4 / pi / const.epsil0 *...
        q(type_i) * q(type_j) / r; % J
end

function F = Force_Coulombic(const, bias, pos_i, pos_j, type_i, type_j)
%     F = 0;
    q = const.q + bias(:,1)';
    
    r = pos_i-pos_j;
    r = sqrt(sum(r.^2)); % A
    r = r * 1e-10; % m
    dz = pos_i(3)-pos_j(3); % A
    dz = dz * 1e-10; % m
    
    F = const.e^2 / 4 / pi / const.epsil0 *...
        q(type_i) * q(type_j) / -r^3 * dz; % N
end

function E = Energy_vdW(const, bias, pos_i, pos_j, type_i, type_j)
%     E = 0;
    D = const.D + bias(:,2)';
    R = const.R + bias(:,3)';
    
    r = pos_i-pos_j;
    r = sqrt(sum(r.^2));
    
    D = sqrt(D(type_i) * D(type_j));
    R = (R(type_i) + R(type_j))/2;
    D = D * const.D_convert; % J
    
    E = D * ( (R/r)^12 - 2*(R/r)^6 ); % J
end

function F = Force_vdW(const, bias, pos_i, pos_j, type_i, type_j)
%     F = 0;
    D = const.D + bias(:,2)';
    R = const.R + bias(:,3)';
    
    r = pos_i-pos_j;
    r = sqrt(sum(r.^2)); % A
    r = r * 1e-10; % m
    dz = pos_i(3)-pos_j(3); % A
    dz = dz * 1e-10; % m
    
    D = sqrt(D(type_i) * D(type_j));
    R = (R(type_i) + R(type_j))/2; % A
    D = D * const.D_convert; % J
    R = R * 1e-10; % m
    
    F = - D * R^12 * 12 / r^14 * dz ...
        - D * R^6 * 6 / r^8 * dz; % N
end

function z = FindZeros(d, mat, init, pos)
    n = size(mat,2);
    z = zeros(1,n);
    for i = 1:n
        c = mat(:,i);
        
%         [xData, yData] = prepareCurveData(d', c');
%         fun = fit(xData, yData, 'splineinterp');
%         
%         z(i) = fzero(fun, init);
% 
        l2 = c(1:end-1) .* c(2:end) < 0;
        ids = find(l2);
        i_init = find(d==init);
        
        if ~isempty(ids)
            dists = abs(ids-i_init);
            iids = find(dists == min(dists),1);
            ids = ids(iids);

            % refine the force zero point, linear interpolation
            assert(c(ids)*c(ids+1)<0)

            slope = (c(ids) - c(ids+1)) / (d(ids) - d(ids+1));
            Fz = d(ids) - c(ids)/ slope;
            z(i) = Fz;

        else
            disp(['Zero point not found!  [', num2str(pos(1)), ',', num2str(pos(2)),']'])
            z(i) = init;
            figure(1)
                scatter(pos(1), pos(2),50,'rx')
            figure(2); hold on;
                plot(c)
%                 pause();
        end
        
    end
end

function [F, Fz] = ForceCalc(pos, const, bias, atom_pos, cutoff_range_A, convert, d_list, init)
    F = zeros(length(d_list),1);
    
    % select atoms in cutoff distance (3d)
    select = atom_pos(:,2:4) - [pos,0]; % A
    select = sqrt(sum(select.^2,2)) < cutoff_range_A;
    pos2 = atom_pos(select,:); % list for neighboring atoms
    pos2(:,2:3) = pos2(:,2:3) - pos; % xy relative position in A

    F_coulomb = @(pos_i, pos_j, type_j) Force_Coulombic(const, bias, pos_i, pos_j, 4, type_j);
    F_vdW = @(pos_i, pos_j, type_j) Force_vdW(const, bias, pos_i, pos_j, 4, type_j);
    F_total = @(pos_i, pos_j, type_j) F_coulomb(pos_i, pos_j, type_j) + F_vdW(pos_i, pos_j, type_j);

    for i = 1:length(d_list)
        z = d_list(i);

            for j = 1:size(pos2,1)
                F(i) = F(i) + F_total(pos2(j,2:4), [0,0,z] + Tip(i_tip,:), pos2(j,1));
            end
    end
    
    Fz = FindZeros(d_list, F, init, pos);
end


