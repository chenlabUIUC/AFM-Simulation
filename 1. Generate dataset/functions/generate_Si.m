%% Function to generate Si atom of different facet and size
% Generate Si atom positions controled by radius, height, and whether it's
% (111) facet or (100) facet.

function atomic_positions = generate_Si(radius, height, lattice_points_per_side, rot111)
    lattice_constant = 5.43;  % Lattice constant of silicon in Angstroms
%     radius = 3;  % Radius of the cylindrical region
%     height = 6;  % Height of the cylindrical region
%     lattice_points_per_side = 6;  % Number of unit cells per side

    atomic_positions = generate_silicon_positions_in_cylinder(radius, height, lattice_constant, lattice_points_per_side, rot111);
end


function atomic_positions = generate_silicon_positions_in_cylinder(radius, height, lattice_constant, lattice_points_per_side, rot111)
    % Define the unit cell coordinates of silicon atoms in fractional coordinates
    fractional_coords = [
        0.00 0.00 0.00;
        0.50 0.50 0.00;
        0.00 0.50 0.50;
        0.50 0.00 0.50;
        0.25 0.25 0.25;
        0.75 0.75 0.25;
        0.75 0.25 0.75;
        0.25 0.75 0.75;
    ];
    
    % Generate all atomic positions in the cubic lattice
    atomic_positions = [];
    for x = -lattice_points_per_side:lattice_points_per_side
        for y = -lattice_points_per_side:lattice_points_per_side
            for z = -lattice_points_per_side:lattice_points_per_side
                fractional_translated = fractional_coords + [x y z];
                atomic_positions = [atomic_positions; fractional_translated];
            end
        end
    end
    
    % Convert fractional coordinates to Cartesian coordinates using the lattice constant
    atomic_positions = atomic_positions * lattice_constant;
    
    if rot111
        rot_z = @(x) [cosd(x), -sind(x), 0;
                      sind(x), cosd(x), 0;
                      0,0,1];
        rot_y = @(x) [cosd(x), 0, -sind(x);
                      0,1,0;
                      sind(x), 0, cosd(x);];

        atomic_positions = atomic_positions * rot_z(45);
        atomic_positions = atomic_positions * rot_y(acos(1/sqrt(3))*180/pi);
    end
    
    % Filter out atoms outside the cylindrical region
    x = atomic_positions(:, 1);
    y = atomic_positions(:, 2);
    z = atomic_positions(:, 3);
    distances = sqrt(x.^2 + y.^2);
    mask = distances <= radius & z <= height & z>=-1e-5;
    atomic_positions = atomic_positions(mask, :);
end


function plot_atoms_with_bonds(atomic_positions, lattice_constant)
    % Number of atoms
    num_atoms = size(atomic_positions, 1);

    % Plot atomic positions as scatter points
    scatter3(atomic_positions(:, 1), atomic_positions(:, 2), atomic_positions(:, 3), 'filled');
    hold on;

    % Plot bonds between neighboring atoms
    for i = 1:num_atoms
        for j = i+1:num_atoms
            distance = norm(atomic_positions(i, :) - atomic_positions(j, :));
            if distance <= 1.001*sqrt(3)/4 * lattice_constant  % You can adjust the bond length threshold as needed
                plot3([atomic_positions(i, 1), atomic_positions(j, 1)], ...
                      [atomic_positions(i, 2), atomic_positions(j, 2)], ...
                      [atomic_positions(i, 3), atomic_positions(j, 3)], 'k-');
            end
        end
    end

    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    axis equal;

    grid off;
    hold off;
end

