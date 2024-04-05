%% AtomicSimulation_defect.m
% Adding defects to the perfect lattices with the following workflow:
% 1. Delete Al atoms (define by bonds)
% 1.1. Add bonds, prepare neighboring list
% 1.2. Define deleted Al atoms
% 1.3. Delete related O-H
% 1.4. Identify edge Al atoms for further add Al atoms
% 2. Add Al atoms (library style)
% 2.1. Make library of adatoms (from DFT, make into library style)
% 2.2. Randomly assign leaving atom along the edge
% 2.3. Substitute perfect structure to adatom library
%
% INPUT:
%   Bond thresholds: 'OH-thresh','AlO_thresh','AlAl_thresh'
%   Pixel size: 'convert', Å/px
%   Randomizing parameters:
%       'Init_vacancy_prob': initially generating the vacancy randomly.
%       'poisson_lambda': assume vacancy size as random poisson
%                         distribution.
%       'large_vac_count': besides small vacancies control by
%                          posson distribution, we generate some bigger
%                          vacancies.
%       'large_vac_size': control the size of the big vacancies.
%       'Leave_zigzag_prob': probability of a zigzag state Al atom to be
%                            replaced.
%       'Leave_armchair_twoAtom_prob': probability of a pair of armchair
%                         state Al atoms to be replaced as two-atom etching
%                         case.
%       'Leave_armchair_oneAtom_prob': probability of a pair of armchair
%                         state Al atoms to be replaced as one-atom etching
%                         case.
%   'plot_bonds': 1*3 logical array, control whether to plot O-H, Al-O and
%                 Al-Al bonds in the final structure
%
% OUTPUT:
% 'output' structure with 2 fields:
%   'output.atom_pos': n*4 matrix, with 1st column as atom type (1 for Al, 2
%                      for O, 3 for H) and 2nd-4th columnns as atom positions
%                      in Å.
%   'output.atom_size': for further plotting
%
% HISTORY:  Written by Chang (Chris) Qian Last modified by Chang (Chris)
% Qian on 04/04/2024


clear; close all; clc;
addpath('functions/')
rng(1)

% For demonstration
num_structure = 1; % How many perfect structures (corresponding to different in-plane rotation)
num_generate = [1]; % How many defect structure to generate for each perfect structure

% % For real training set generation
% num_structure = 7;
% num_generate = [5,10,20,40,20,10,5]; % 7 structures


for structure_ind = 1:num_structure
for repeating = 1:num_generate(structure_ind)
load(['output/Perfect structure/Perfect structures_',num2str(structure_ind),'.mat'])
atom_pos = output.atom_pos;
atom_size = output.atom_size;

%%%% workflow control %%%%
plot_in_A = 1;

save_output = 0;
save_name = ['output/Defected structure/Defected structures_',...
             num2str(structure_ind),'-',num2str(repeating),'.mat'];

OH_thresh = 1.3; % Å
AlO_thresh = 2.5; % Å
AlAl_thresh = 3.3; % Å

plot_size_px = 200;
convert = 200/256; % Å/px

Init_vacancy_prob = 0.005 + rand()*0.015; % 0.5-2.0%
poisson_lambda = 2 + rand(); % 2-3
large_vac_count = randi([8,15]); % 8-15
large_vac_size = [20,150];

Leave_zigzag_prob = 0.4;
Leave_armchair_twoAtom_prob = 0.7;
Leave_armchair_oneAtom_prob = 0.8; % if not two atom etching.

plot_bonds = [1,1,1];
skip_plotting = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%

% load adatom library
load('Adatom Library_20240313.mat');


%% 1. Delete Al atoms (define by bonds)
%% 1.1. Add bonds, neighboring list
[AlAl_bonds, AlO_bonds, OH_bonds] = makeBonds(atom_pos, AlAl_thresh, AlO_thresh, OH_thresh);


%% 1.2. Define deleted Al atoms

% create initial seed
Al_inds = find(atom_pos(:,1)==1);
Inds_inv = zeros(size(atom_pos,1),1);
    for i = 1:size(Al_inds,1)
        Inds_inv(Al_inds(i)) = i;
    end
Al_delete = rand(size(Al_inds)) < Init_vacancy_prob;

inds_temp = find(Al_delete);
sel = atom_pos(Al_inds(Al_delete),4) < -2.5; % 2nd layer
Al_delete(inds_temp(sel)) = 0;

% Al neighbor list
Al_neighbors = zeros(size(Al_inds,1), 3);
for i = 1:size(Al_inds,1)
    ind = Al_inds(i);
    sel1 = AlAl_bonds(AlAl_bonds(:,1)==ind,2);
    sel2 = AlAl_bonds(AlAl_bonds(:,2)==ind,1);
    sel = [sel1',sel2'];
    Al_neighbors(i,1:length(sel)) = sel;
end

% define area size
Al_vacancy_size = poissrnd(poisson_lambda,size(Al_inds,1),1);
Al_vacancy_size(~Al_delete) = 0;
Al_delete_site_inds = find(Al_delete);
randomIndices = Al_delete_site_inds(randperm(length(Al_delete_site_inds), large_vac_count));
randomSizes = randi(large_vac_size, large_vac_count,1);
Al_vacancy_size(randomIndices) = randomSizes;
Al_delete_site_inds = cat(2, Al_delete_site_inds, Al_vacancy_size(Al_delete_site_inds));

% expand the vacancy
for i = 1:size(Al_delete_site_inds,1)
    del_list = [Al_delete_site_inds(i,1)];
    del_num = Al_delete_site_inds(i,2);
    
    if del_num<2; continue; end
    
    while del_num>0 && ~isempty(del_list)
        temp = del_list(1);
        del_list(1) = [];
        del_num = del_num - 1;
        Al_delete(temp) = 1;
        
        % append neighbors, that are not deleted.
        n_temp = Al_neighbors(temp,:);
        n_temp = n_temp(n_temp~=0);
        n_temp = Inds_inv(n_temp);
        n_temp = n_temp(~Al_delete(n_temp));
        n_temp = reshape(n_temp, 1, []);
        del_list = [del_list, n_temp];
    end
end

% fill the single vacancies back, because this won't delete any O-H, not
% physical.
Al_delete_site_inds = find(Al_delete);
for i = 1:size(Al_delete_site_inds,1)
    temp = Al_delete_site_inds(i);
    temp = Al_neighbors(temp,:);
    temp = temp(temp~=0);
    temp = Al_delete(Inds_inv(temp));
    if sum(temp)==0
        Al_delete(Al_delete_site_inds(i)) = 0;
    end
end
disp(['Total deleted Al atom: ',num2str(sum(Al_delete))]);

%% 1.3. Identify edge Al atoms for further add Al atoms.
% Delete dangling bonds?
delete_dangling = 1; % number of iteration to delete dangling

Al_edge_type = zeros(size(Al_inds));
% 0: in bulk; >0: at/near edge; -1:deleted
% 1: zigzag
% 2: armchair
% 3: near edge (anchors)
% 4: dangling
% 5: bridging
for i = 1:size(Al_inds,1) % Determine -1, 4. 0 as bulk, 1 as defected.
    if Al_delete(i); Al_edge_type(i) = -1; end
    nn_temp = Al_neighbors(i,:);
    nn_temp = nn_temp(nn_temp~=0);
    nn_temp = Inds_inv(nn_temp);
    n_deleted_temp = sum(Al_delete(nn_temp)==0);
    
    switch n_deleted_temp
        case 1 % dangling
            Al_edge_type(i) = 4;
        case 2 % zigzag, armchair, bridging. Wait for categorization.
            Al_edge_type(i) = 1;
        case 3 % bulk
            continue
    end
end

if delete_dangling
    while delete_dangling>0
        delete_dangling = delete_dangling-1;
        Al_delete(Al_edge_type==4) = 1;

        Al_edge_type = zeros(size(Al_inds));
        for i = 1:size(Al_inds,1) % Determine -1, 4. 0 as bulk, 1 as defected.
            if Al_delete(i); Al_edge_type(i) = -1; continue; end
            nn_temp = Al_neighbors(i,:);
            nn_temp = nn_temp(nn_temp~=0);
            nn_temp = Inds_inv(nn_temp);
            n_deleted_temp = sum(Al_delete(nn_temp)==0);

            switch n_deleted_temp
                case 1 % dangling
                    Al_edge_type(i) = 4;
                case 2 % zigzag, armchair, bridging. Wait for categorization.
                    Al_edge_type(i) = 1;
                case 3 % bulk
                    continue
            end
        end
    end
end

list_temp = find(Al_edge_type==0);
for ind = 1:size(list_temp,1) % Determine 0, 3
    i = list_temp(ind);
    nn_temp = Al_neighbors(i,:);
    nn_temp = nn_temp(nn_temp~=0);
    nn_temp = Inds_inv(nn_temp);
    nn_temp = Al_edge_type(nn_temp);
    nn_temp = nn_temp(nn_temp~=-1);
    if sum(nn_temp==1 | nn_temp==4)>0; Al_edge_type(i) = 3; end % near edge.
end

list_temp = find(Al_edge_type==1);
for ind = 1:size(list_temp,1) % Determine 1, 5
    i = list_temp(ind);
    nn_temp = Al_neighbors(i,:);
    nn_temp = nn_temp(nn_temp~=0);
    nn_temp = Inds_inv(nn_temp);
    nn_temp = Al_edge_type(nn_temp);
    nn_temp = nn_temp(nn_temp~=-1);
    
    if sum(nn_temp==3)==2; Al_edge_type(i) = 1; end % zigzag case
    if sum(nn_temp==3)==1; Al_edge_type(i) = 2; end % armchair and bridge
    if sum(nn_temp==3)==0; Al_edge_type(i) = 5; end % bridge only (>2 atom bridge)
end

list_temp = find(Al_edge_type==2);
for ind = 1:size(list_temp,1) % Determine 2, 5
    i = list_temp(ind);
    nn_temp = Al_neighbors(i,:);
    nn_temp = nn_temp(nn_temp~=0);
    nn_temp = Inds_inv(nn_temp);
    nn_temp = Al_edge_type(nn_temp);
    nn_temp = nn_temp(nn_temp~=-1);
    
    if sum(nn_temp==5)>0; Al_edge_type(i) = 5; end % bridge only (>2 atom bridge)
    
    if sum(nn_temp==2)>0; Al_edge_type(i) = 2; end % armchair and bridge (2 atom bridge)
end

%% 1.4. Delete related Al-O-H
O_inds = find(atom_pos(:,1)==2);
    for i = 1:size(O_inds,1)
        Inds_inv(O_inds(i)) = i;
    end % Inds_inv
O_delete = O_inds==0; % all logical 0
H_inds = find(atom_pos(:,1)==3);
    for i = 1:size(H_inds,1)
        Inds_inv(H_inds(i)) = i;
    end % Inds_inv
H_delete = H_inds==0; % all logical 0
% Check oxygen atoms
for i = 1:size(O_inds,1)
    i_temp = O_inds(i);
    n_temp = AlO_bonds(AlO_bonds(:,2)==i_temp,1);
    n_temp = Inds_inv(n_temp);
    if sum(Al_delete(n_temp)==0)==0 % all connected Al atom deleted
        O_delete(i) = 1;
    end
end
% Check hydrogen atoms
for i = 1:size(H_inds,1)
    i_temp = H_inds(i);
    n_temp = OH_bonds(OH_bonds(:,2)==i_temp,1);
    n_temp = Inds_inv(n_temp);
    if sum(O_delete(n_temp)==0)==0 % all connected O atom deleted
        H_delete(i) = 1;
    end
end

% Summarize all deleted atoms
atoms_delete = zeros(size(Inds_inv));
atoms_delete(Al_inds(Al_delete)) = 1;
atoms_delete(O_inds(O_delete)) = 1;
atoms_delete(H_inds(H_delete)) = 1;

% plotAtomPos(atom_pos(~atoms_delete,:),4, plot_size_px, convert, atom_size)


%% 2. Add Al atoms (library style)
%% 2.1. Make library of adatoms (from DFT, make into library style)
% 1-8: armchair two atom
% 9-13: armchair one atom
% 14-18: zigzag one atom

%% 2.2. Randomly assign leaving atom along the edge
% Atom groups to leave: 'leaving_groups'
ind_defect = Al_inds(Al_edge_type==1 | Al_edge_type==2);
num_defect = length(ind_defect);
leaving_groups = cell(num_defect,1);

for idx = 1:size(ind_defect,1) % Preparing leaving_groups
    i = ind_defect(idx);
    atoms_in_group = [i];
    % Al
    inds_temp = Al_neighbors(Inds_inv(i),:);
    inds_temp = inds_temp(inds_temp~=0);
    inds_temp = inds_temp(Al_delete(Inds_inv(inds_temp))==0);
    atoms_in_group = [atoms_in_group, inds_temp];
    % O
    for ind = atoms_in_group
        O_temp = AlO_bonds(AlO_bonds(:,1)==ind,2);
        for O_ind_temp = O_temp'
            Al_temp = AlO_bonds(AlO_bonds(:,2)==O_ind_temp,1);
            Al_temp = Al_temp(~atoms_delete(Al_temp));
            if sum(ismember(Al_temp, atoms_in_group)==0)==0
                atoms_in_group = [atoms_in_group, O_ind_temp];
            end
        end
    end
    % H
    for ind = atoms_in_group(4:end)
        H_temp = OH_bonds(OH_bonds(:,1)==ind,2);
        atoms_in_group = [atoms_in_group, H_temp];
    end
    leaving_groups{idx} = atoms_in_group;
end 

% Index in library to replace: 'replacing_ind'
replacing_ind = zeros(num_defect,1);
nn_armchair_sites = zeros(num_defect,1);
for i = 1:size(ind_defect,1) % Preparing replacing_ind
    defect_type_temp = Inds_inv(ind_defect(i));
    defect_type_temp = Al_edge_type(defect_type_temp);
    
    if defect_type_temp == 1 % zigzag
        if rand() < Leave_zigzag_prob
            replacing_ind(i) = randi([14,18]);
        end
    elseif defect_type_temp == 2 % armchair
        % site not replaced
        sel1 = replacing_ind(i)==0; 
        
        nn_armchair_site = Al_neighbors(Inds_inv(ind_defect(i)),:);
        nn_armchair_site = nn_armchair_site(nn_armchair_site~=0);
        nn_armchair_site = Inds_inv(nn_armchair_site);
        nn_armchair_site = nn_armchair_site(Al_edge_type(nn_armchair_site)==2);
        
        if isempty(nn_armchair_site)
            continue; % it's a bridge near dangling Al atom
        else
            % neighboring armchair site not replaced
            nn_armchair_sites(i) = find(ind_defect==Al_inds(nn_armchair_site));
            sel2 = replacing_ind(ind_defect==Al_inds(nn_armchair_site))==0; 
        end
        
        if ~sel1; continue; end
        
        if rand() < Leave_armchair_twoAtom_prob && sel2 % Two atom
            random_ind = randi([1,8]);
            replacing_ind(i) = random_ind;
            replacing_ind(ind_defect==Al_inds(nn_armchair_site)) = random_ind;
        elseif rand() < Leave_armchair_oneAtom_prob % One atom
            replacing_ind(i) = randi([9,13]);
        end
    end
end

%% 2.3. Substitute perfect structure to adatom library
% replace atoms except 2 anchor Al atoms.
Add_atoms = zeros(0,4);
% replacing_ind(nn_armchair_sites>0 & nn_armchair_sites<[1:size(nn_armchair_sites,1)]) = 0;

for i = 1:size(replacing_ind,1)
    if replacing_ind(i) == 0; continue; end
    
    adding_atoms = adatom_lib{replacing_ind(i)}.atom_pos;
    
    if replacing_ind(i) > 8 % zigzag one atom, or armchair one atom
        anchor = leaving_groups{i}(2:3);
        leaving = leaving_groups{i}(1);        
    else % armchair two atom
        anchor = leaving_groups{i}(2:3);
        leaving = leaving_groups{i}(1);
        nn_armchair_site = nn_armchair_sites(i);
        if nn_armchair_site > i
            continue; 
        end % already replaced.
        
        anchor = [anchor, leaving_groups{nn_armchair_site}(2:3)];
        leaving = [leaving, leaving_groups{nn_armchair_site}(1)];
        sel = ismember(anchor, leaving);
        anchor(sel)=[];
    end
    
    v = atom_pos(anchor,2:3);
    v = v - atom_pos(leaving(1),2:3);
    if det(v)<0; anchor = flip(anchor); end

    sourcePoints = [adding_atoms(1:2,2:3)];
    targetPoints = [atom_pos(anchor,2:3)];

    tform = fitgeotrans(sourcePoints, targetPoints, 'nonreflectivesimilarity');
    scale = det(tform.T(1:2,1:2));
    if scale>1.2
        % disp(['Bridging case, no DFT structure available'])
        continue
    end
    adding_atoms(:,2:3) = transformPointsForward(tform, adding_atoms(:,2:3));

    atoms_delete(leaving_groups{i}([1,4:end])) = 1;
    if length(leaving)>1
        atoms_delete(leaving_groups{nn_armchair_site}([1,4:end])) = 1;
    end
    Add_atoms = cat(1, Add_atoms, adding_atoms(3:end,:));
    
end

%% Summarize, delete and adding, redo bond network.
atom_pos(logical(atoms_delete),:) = [];
atom_pos = cat(1, atom_pos, Add_atoms);

if ~skip_plotting 
[AlAl_bonds, AlO_bonds, OH_bonds] = makeBonds(atom_pos, AlAl_thresh, AlO_thresh, OH_thresh);

% plotAtomPos(atom_pos, 5, plot_size_px, convert, atom_size)
%% Draw crystal
figure(1); clf; hold on; axis equal;

if plot_in_A
    atom_pos_plot = atom_pos;
    axis_lim = plot_size_px * convert;
else
    atom_pos_plot = atom_pos;
    atom_pos_plot(:,2:4) = atom_pos_plot(:,2:4)./convert;
    axis_lim = plot_size_px;
end

axis_lim = axis_lim;
% sel = atom_pos_plot(:,1)==1 &...
%       atom_pos_plot(:,4) > -2.5 & ...
%       atoms_delete==0;
% sel = atom_pos_plot(:,1)==1 &...
%       atom_pos_plot(:,4) > -2.5;
sel = atom_pos_plot(:,2)>0 & atom_pos_plot(:,2)<axis_lim &...
      atom_pos_plot(:,3)>0 & atom_pos_plot(:,3)<axis_lim &...
      atom_pos_plot(:,4) > -2.5;
  
% Atoms
scatter3(atom_pos_plot(sel,2),atom_pos_plot(sel,3),atom_pos_plot(sel,4),...
         atom_size(atom_pos(sel,1)).^2*30, atom_pos_plot(sel,4),'filled')
colormap('cool')

% % Bonds
if plot_bonds(1) % Al-Al
for i = 1:size(AlAl_bonds,1)
    if sum(sel(AlAl_bonds(i,:))==0)==0
    plot3(atom_pos(AlAl_bonds(i,:),2),...
          atom_pos(AlAl_bonds(i,:),3),...
          atom_pos(AlAl_bonds(i,:),4),'k-')
    end
end
end 
if plot_bonds(2) % Al-O
for i = 1:size(AlO_bonds,1)
    if sum(sel(AlO_bonds(i,:))==0)==0
    plot3(atom_pos_plot(AlO_bonds(i,:),2),...
          atom_pos_plot(AlO_bonds(i,:),3),...
          atom_pos_plot(AlO_bonds(i,:),4),'k-')
    end
end
end 
if plot_bonds(3) % O-H
for i = 1:size(OH_bonds,1)
    if sum(sel(OH_bonds(i,:))==0)==0
    plot3(atom_pos_plot(OH_bonds(i,:),2),...
          atom_pos_plot(OH_bonds(i,:),3),...
          atom_pos_plot(OH_bonds(i,:),4),'k-')
    end
end
end

axis equal
% caxis([-1, 1] .* 0.05)
xlim([0,axis_lim])
ylim([0,axis_lim])


% defect type
for i = 0:-1
    logicalArray = false(size(atom_pos,1),1);
    logicalArray(Al_inds(Al_edge_type==i)) = true;
    switch i
        case 0
%             scatter(atom_pos(sel & logicalArray,2),...
%                     atom_pos(sel & logicalArray,3),...
%                     'ko','filled');
        case 1
            scatter3(atom_pos(sel & logicalArray,2),...
                    atom_pos(sel & logicalArray,3),...
                    atom_pos(sel & logicalArray,4),...
                    'ro','filled');
        case 2
            scatter3(atom_pos(sel & logicalArray,2),...
                    atom_pos(sel & logicalArray,3),...
                    atom_pos(sel & logicalArray,4),...
                    'bo','filled');
        case 3
            scatter3(atom_pos(sel & logicalArray,2),...
                    atom_pos(sel & logicalArray,3),...
                    atom_pos(sel & logicalArray,4),...
                    'go','filled');
        case 4
            scatter3(atom_pos(sel & logicalArray,2),...
                    atom_pos(sel & logicalArray,3),...
                    atom_pos(sel & logicalArray,4),...
                    'co','filled');
        case 5
            scatter3(atom_pos(sel & logicalArray,2),...
                    atom_pos(sel & logicalArray,3),...
                    atom_pos(sel & logicalArray,4),...
                    'mo','filled');
    end
end

end

%% Prepare output file
output = struct;
output.atom_pos = atom_pos;
output.atom_size = atom_size;
if save_output
    save(save_name, 'output');
end

end
end

%% Functions
function [AlAl_bonds, AlO_bonds, OH_bonds] = makeBonds(atom_pos, AlAl_thresh, AlO_thresh, OH_thresh)
    AlAl_bonds = [];
    AlO_bonds = [];
    OH_bonds = [];
    
    t3 = size(atom_pos,1);
    for i = 1:size(atom_pos,1)
        if i==1; tic; end
        
        d = atom_pos(i+1:end, 2:4) - atom_pos(i, 2:4);
        d = sqrt(sum(d.^2, 2));
        id = atom_pos(i+1:end, 1);
        ind = [i+1:size(atom_pos,1)];

        switch atom_pos(i,1)
            case 1
                sel = id==1 & d<AlAl_thresh;
                inds = ind(sel);
                inds2 = ones(size(inds)) * i;
                AlAl_bonds = [AlAl_bonds; inds2', inds'];

                sel = id==2 & d<AlO_thresh;
                inds = ind(sel);
                inds2 = ones(size(inds)) * i;
                AlO_bonds = [AlO_bonds; inds2', inds'];

            case 2
                sel = id==1 & d<AlO_thresh;
                inds = ind(sel);
                inds2 = ones(size(inds)) * i;
                AlO_bonds = [AlO_bonds; inds', inds2'];

                sel = id==3 & d<OH_thresh;
                inds = ind(sel);
                inds2 = ones(size(inds)) * i;
                OH_bonds = [OH_bonds; inds2', inds'];

            case 3
                sel = id==2 & d<OH_thresh;
                inds = ind(sel);
                inds2 = ones(size(inds)) * i;
                OH_bonds = [OH_bonds; inds', inds2'];

            if i==1
                t = toc;
                t2 = t * t3;
                disp([num2str(t,'%.3e'), 's  ', 'Estimated: ',...
                      num2str(t2, '%.3e'), 's'])
            end
        end
    end
end

function plotAtomPos(atom_pos, ii, plot_size_px, convert, atom_size)
    figure(ii); clf; hold on; axis equal;

    atom_pos_plot = atom_pos;
    axis_lim = plot_size_px * convert;
    
    sel = atom_pos_plot(:,2)>0 & atom_pos_plot(:,2)<axis_lim &...
          atom_pos_plot(:,3)>0 & atom_pos_plot(:,3)<axis_lim &...
          atom_pos_plot(:,4) > -2.5;

    % Atoms
    scatter3(atom_pos_plot(sel,2),atom_pos_plot(sel,3),atom_pos_plot(sel,4),...
             atom_size(atom_pos(sel,1)).^2*30, atom_pos_plot(sel,4),'filled')
    colormap('cool')
end

