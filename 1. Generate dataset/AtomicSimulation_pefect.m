%% AtomicSimulation_perfect.m
% Read .vasp file of gibbsite. Tile the crystal over the whole image
% (256*256 px), with extended buffer area (50 Å).
% 
% INPUT:
% 'poscar_file': .vasp file of crystal structure
% 'plot_size_px': size of the AFM image in pixel
% 'plot_buffer_A': size of buffer area in Å
% 'convert': pixel size of the AFM image
% 
% OUTPUT:
% 'output' structure with 3 fields:
%   'output.atom_pos': n*4 matrix, with 1st column as atom type (1 for Al, 2
%                      for O, 3 for H) and 2nd-4th columnns as atom positions
%                      in Å.
%   'output.atom_size': for further plotting
%   'output.atom_color': for further plotting
% 
% HISTORY:  Written by Chang (Chris) Qian Last modified by Chang (Chris)
% Qian on 04/04/2024

clear; close all; clc;
addpath('functions/')

% rot_offset_list = [-3:3]/2; % For real dataset, we need to sample lattice of different orientation.
rot_offset_list = [0]; % In this example, we just generate one.

for rot_offset_ind = 1:length(rot_offset_list)

rot_offset = rot_offset_list(rot_offset_ind);
    
%%%% workflow control %%%%
plot_in_A = 1;

save_output = 1;
save_name = ['output/Perfect structure/Perfect structures_',...
             num2str(rot_offset_ind),'.mat'];

poscar_file = 'gibbsite.vasp';
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_size_px = 256; % px
plot_buffer_A = 50; % Å

z_range = [-0.7,0.3]; % This range is for 2 layer of gibbsite structure.

rot = 90 + 80.243 + rot_offset; % 2d rotation, deg.
convert = 200/256; % Pixel size, Å/px
%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_size_A = plot_size_px * convert;
plot_buffer_px = plot_buffer_A / convert; % px

%% Constants
rot_mat = RotMat_handle();
Const = PrepConst;

%% Parameters in POSCAR
f = import_poscar(poscar_file);

% Atom, read in
pos = f.coords;
symbols = f.symbols;
atom_count = f.atomcount;

% Box, read in
lattice = f.lattice;

% Id, size list
assert(size(pos,1) == sum(atom_count));
atom_id = [];
atom_size = [];

for i = 1:length(symbols)
    atom_id = cat(1,atom_id,ones(atom_count(i),1)*i);
end
for i = 1:length(symbols)
    atom_size(i)=atomic_radius(symbols{i});
end

atom_color = lines(length(symbols));

%% Tile the lattice over given viewing area, plus buffer area
% Box, shifted
L = plot_size_px + plot_buffer_px*2; % in px
L = L * convert; % in Å
extrema = [L,0; L,L; 0,L; 0,0];

% Rotated lattice
a1_rotated = lattice(1,:) * rot_mat(rot);
a2_rotated = lattice(2,:) * rot_mat(rot);

% Tiling repetition, in 2d
N1 = extrema * a1_rotated(1:2)' ./ norm(a1_rotated)^2;
N2 = extrema * a2_rotated(1:2)' ./ norm(a2_rotated)^2;
N1 = [floor(min(N1)), ceil(max(N1))]; n1 = N1(2)-N1(1)+1;
N2 = [floor(min(N2)), ceil(max(N2))]; n2 = N2(2)-N2(1)+1;

% Repetition in z
N3 = floor(z_range(1)):ceil(z_range(2))-1; n3 = max(N3)-min(N3)+1;

% Generate all atomic positions
N_atom_cell = size(pos,1);
atom_pos = zeros(n1*n2*n3*N_atom_cell, 4);

counter = 0;
for i = N1(1):N1(2)
for j = N2(1):N2(2)
for k = N3
    pos2 = (pos + [i,j,k]) .* [lattice(1,1),lattice(2,2),lattice(3,3)];
    pos2 = cat(2, atom_id, pos2);
    atom_pos(counter*N_atom_cell+[1:N_atom_cell], :) = pos2;
    counter = counter+1;
end; end; end

% Rotate
atom_pos(:,2:4) = atom_pos(:,2:4) * rot_mat(rot);

% Trim in x,y,z
select1 = atom_pos(:,2)<0 | atom_pos(:,2)>L;
select2 = atom_pos(:,3)<0 | atom_pos(:,3)>L;
select3 = atom_pos(:,4)<z_range(1)*lattice(3,3) | atom_pos(:,4)>z_range(2)*lattice(3,3);
atom_pos(select1 | select2 | select3,:) = [];

% Shift the buffer out
atom_pos(:,2:3) = atom_pos(:,2:3) - [plot_buffer_px, plot_buffer_px].*convert;

%% Draw crystal
figure(1); clf; hold on; axis equal;

if plot_in_A
    scatter3(atom_pos(:,2),atom_pos(:,3),atom_pos(:,4),...
             atom_size(atom_pos(:,1)).^2*30, atom_pos(:,4),'filled')
    colormap('cool')
    axis equal

    caxis([-1, 1] .* 0.05)
    xlim([0,plot_size_px].*convert)
    ylim([0,plot_size_px].*convert)
else
    scatter3(atom_pos(:,2)./convert,atom_pos(:,3)./convert,atom_pos(:,4)./convert,...
             atom_size(atom_pos(:,1)).^2*30, atom_pos(:,4),'filled')
    colormap('cool')
    axis equal

    caxis([-1, 1] .* 0.05)
    xlim([0,plot_size_px])
    ylim([0,plot_size_px])
end

% box
% plotBox(lattice, 'k-');

%% Prepare output file
output = struct;
output.atom_pos = atom_pos;
output.atom_size = atom_size;
output.atom_color = atom_color;
if save_output
    save(save_name, 'output');
end

end

%% Functions
function plotBox(lattice, varargin)
    a = lattice(1,:); b = lattice(2,:); c = lattice(3,:);
    coord = [0,0,0; a; a+b; b;
             0,0,0; c; c+a; c+a+b; c+b; c;
             0,0,0; b; b+c; a+b+c; a+b; a; a+c];
	plot3(coord(:,1),coord(:,2),coord(:,3), varargin{:});
end
