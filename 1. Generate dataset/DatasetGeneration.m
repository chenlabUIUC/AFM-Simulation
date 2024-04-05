%% DatasetGeneration_demo.m
% From force matrix, generate dataset including training images and ground
% truth. 
%
% INPUT:
%   Number of images to generate:
%       'n_train' and 'n_val' for training and validation set.
%   Random parameters related to the image generation:
%       'params' as a structure, define from line 114, including parameters
%       for Si atom interaction, tip size and shape, random noise level,
%       intensity clipping, z-map adjusting, random intensity patches.
%   Parameters related to ground truth (GT) generation:
%       'Params.unmoved_Al': z height threshold to determine if the Al atom
%       is original Al atom or adatom.
%       'Params.Al_intensity': intensity of the Gaussian peak of Al atom in
%       GT.
%       'Params.Al_sigma': sigma of the Gaussian peak of Al atom in GT.
% 
% OUTPUT:
%   'image_stack.h5': training and validation dataset with GT for U-Net
%   training.
%   'params_rec.mat': record of the random parameters for backing out image
%   parameters if necessary.
% 
% HISTORY:  Written by Chang (Chris) Qian Last modified by Chang (Chris)
% Qian on 04/04/2024

clear; close all; clc;
addpath('functions/')

rng(1)
global F_mat atom_pos Params Const

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset_dir = 'output/Dataset/';
force_mat_dir = 'output/Force matrix/';
atom_pos_dir = 'output/Defected structure/';
n_train = 500;
n_val = 20;

Params = struct;
Params.generation_size = 128; % px
Params.convert = 200/256; % Å/px

% workflow control
Params.scatter_check = 1;
Params.check_GT = 1;
Params.check_generated = 1;
Params.verbal = 20;

% More or less constant
Params.unmoved_Al = 0.025;
Params.Al_intensity = 1;
Params.Al_sigma = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

force_mat_list = dir(force_mat_dir);
force_mat_list([force_mat_list(:).isdir]) = [];
atom_pos_list = dir(atom_pos_dir);
atom_pos_list([atom_pos_list(:).isdir]) = [];
n_files = length(force_mat_list);

% number of item generation per file
inds_temp = randi([1,n_files],n_train + n_val,1);
isval = zeros(n_train + n_val,1);
isval(randperm(n_train + n_val,n_val)) = 1;

num_generation = zeros(n_files,1);
for i = 1:n_files
    num_generation(i) = sum(inds_temp == i);
end

% Other initialization
file_counter = 0;
train_counter = 0;
val_counter = 0;
Const = PrepConst();
Const.rot_mat = @(x) [cosd(x), -sind(x), 0;
                      sind(x), cosd(x),  0;
                      0,       0,        1;];
                  
params_rec = [];
image_stack = [];

if Params.check_generated
    figure(1);
    set(gcf,'Position',[100,100,Params.generation_size,Params.generation_size]);
    set(gcf,'Position',[100,100,500,500]);
    set(gca,'Position',[0 0 1 1]);
    figure(2);
    set(gcf,'Position',[150+Params.generation_size,100,Params.generation_size,Params.generation_size]);
    set(gcf,'Position',[700,100,500,500]);
    set(gca,'Position',[0 0 1 1]);
end
                  
if ~isfolder(dataset_dir); mkdir(dataset_dir); end
if ~isfolder([dataset_dir,'img/']); mkdir([dataset_dir,'img/']); end
                  
%% Start generation
for file_ind = 1:n_files
    if num_generation(file_ind)==0; continue; end
    
    % Loading data
    load([atom_pos_dir,atom_pos_list(file_ind).name]);
    atom_pos = output.atom_pos;
    load([force_mat_dir,force_mat_list(file_ind).name]);
    
    for generation_ind = 1:num_generation(file_ind)
        params = struct;
        
        % Si atom
        params.qi = uniformRandomNumber(0.2,2.1); % 2.100
        params.logDi = uniformRandomNumber(-6,-2.5); % -5.7351
        params.Ri = uniformRandomNumber(3,6); % 3.7064
        % Tip atoms/interactions
        params.Rot = uniformRandomNumber(0,360);
        params.Tip_radius = normalRandomNumber(3,0.5);
        params.Tip_height = normalRandomNumber(7,0.5);
        params.SphericalTip = uniformRandomNumber(0.01,0.15);
        params.Decay = uniformRandomNumber(0.1,1.5);
        params.Facet111 = 1;
        params.Tip_atoms = GenerateTipAtoms(params);
        % zero point search
        params.init = 10; % A
        params.d_list = d_list;
        % AFM noise
    	params.Sigma1 = uniformRandomNumber(0.02,0.15);
    	params.drift_prob = [uniformRandomNumber(0.2,0.5),uniformRandomNumber(0.2,0.5)];
    	params.drift_dist = [uniformRandomNumber(1,3),normalRandomNumber(0.4,0.1)]; % px
    	params.Gauss_win = uniformRandomNumber(8,20);
    	params.Gauss_sig = normalRandomNumber(0.5,0.1) * params.Gauss_win;
    	params.Gauss2d = uniformRandomNumber(0.8,1.3);
    	params.Sigma2 = uniformRandomNumber(0.02,0.8); % *sigma_1
        % intensity clip
        params.Autoint = uniformRandomNumber(0.6,2); % A
        params.IntShift = normalRandomNumber(0.5,0.05); % A
        
        % z-map adjust:
        params.adatom_offset = uniformRandomNumber(0.1,1);
        params.vacancy_offset = -uniformRandomNumber(0.1,1);
        params.adatom_sigma = uniformRandomNumber(0,0.5);
        params.vacancy_sigma = uniformRandomNumber(0,0.5);
        
        % random contrast/brightness patch
        params.numPatches = randi([0,8]); % Adjust as needed
        params.maxPatchRadius = 60; % Maximum radius for oval patches
        params.contrastFactorRange = [5/6, 6/5]; % Range of contrast enhancement factors
        params.brightnessFactorRange = [-1.0, 1.0]; % Range of contrast enhancement factors
        params.blendingRadius = 5; % Blending region around patches
        params.patchParams = zeros(params.numPatches, 6);
        cols = 160; rows = 160;
        for i = 1:params.numPatches
            params.patchParams(i,1) = randi(cols);
            params.patchParams(i,2) = randi(rows);
            params.patchParams(i,3) = randi(params.maxPatchRadius);
            params.patchParams(i,4) = randi(params.maxPatchRadius);
            params.patchParams(i,5) = uniformRandomNumber(params.contrastFactorRange(1),params.contrastFactorRange(2));
            params.patchParams(i,6) = uniformRandomNumber(params.brightnessFactorRange(1),params.brightnessFactorRange(2));
        end
        
        % for checking GT
%     	params.Sigma1 = 0.01;
%     	params.Sigma2 = 0.01;
%     	params.Gauss2d = 0.1;
%         params.drift_prob = [0,0];
%         params.Gauss_sig = params.Gauss_win;
        
        %% Generate image
        [I, cm, cM, mat_sz] = ImageGeneration(params);
%         I = flipud(I);
        params.offset_I = [randi([0,size(I,1)-Params.generation_size]),...
                           randi([0,size(I,2)-Params.generation_size])];
%         params.offset_I = [0,0];
        I = I(1+params.offset_I(1):Params.generation_size+params.offset_I(1), ...
              1+params.offset_I(2):Params.generation_size+params.offset_I(2));
        
        if Params.check_generated
            figure(1);
            imshow(I);
            caxis([cm, cM])
        end
        
        %% Generate ground truth (GT)
        offset = [-floor(mat_sz(1)/2)+0.5, -floor(mat_sz(2)/2)+0.5] -...
                  flip(params.offset_I);
        GT = GTGeneration(offset);
%         GT = flipud(GT);
        
        if Params.check_generated
            figure(2);
            imshow(GT);
        end
        
        %% Output
        
        % Match generated image with atomic position
%         figure
%         imshow(flipud(I))
%         hold on
%         sel0 = atom_pos(:,1)==1 & atom_pos(:,4)>Params.unmoved_Al;
%         sel1 = atom_pos(:,1)==1 & abs(atom_pos(:,4))<Params.unmoved_Al;
%         sel = sel0;
%         scatter(atom_pos(sel,2)/Params.convert-floor(mat_sz(1)/2)+0.5,...
%                 atom_pos(sel,3)/Params.convert-floor(mat_sz(2)/2)+0.5,'rx')
%         sel = sel1;
%         scatter(atom_pos(sel,2)/Params.convert-floor(mat_sz(1)/2)+0.5,...
%                 atom_pos(sel,3)/Params.convert-floor(mat_sz(2)/2)+0.5,'b.')
%         caxis([cm,cM])

        GT(:,:,3) = (I-cm) ./ (cM-cm);
        if Params.check_GT
            figure
            imshow(GT)
        end
        
        params_rec = [params_rec, params];
        file_counter = file_counter + 1;
        
        if ~mod(file_counter, Params.verbal)
            disp(['Generated: ',num2str(file_counter)])
        end
        
        if file_counter == 1
            image_stack = GT;
        else
            image_stack = cat(4,image_stack, GT);
        end
        
        imwrite(flipud(GT(:,:,3)), [dataset_dir,'img/',num2str(file_counter),'.tif'])
        
    end
end

image_stack(isnan(image_stack)) = 0;

save([dataset_dir,'params_rec.mat'],'params_rec')
if exist([dataset_dir,'image_stack.h5'],'file')
    delete([dataset_dir,'image_stack.h5']); end
h5create([dataset_dir,'image_stack.h5'], '/data', size(image_stack));
h5write([dataset_dir,'image_stack.h5'], '/data', image_stack);
    
    
function [I, cm, cM, mat_sz] = ImageGeneration(params)
    global F_mat Params Const
    
    %% Image simulation
    params.Tip = [params.qi,...
               10^params.logDi * Const.D_convert,...
               params.Ri];

    % Coulomb
    Force_mat = params.Tip(1) .* F_mat(:,:,:,1);
    % vdW 12
    for i = 1:13
        Force_mat = Force_mat + F_mat(:,:,:,1+i) .* sqrt(params.Tip(2)) .* params.Tip(3)^(13-i);
    end
    % vdW 6
    for i = 1:7
        Force_mat = Force_mat + F_mat(:,:,:,14+i) .* sqrt(params.Tip(2)) .* params.Tip(3)^(7-i);
    end

    Ta = params.Tip_atoms * Const.rot_mat(params.Rot);
        d = sqrt(sum(Ta(:,:).^2,2));
    Ta(:,1:2) = Ta(:,1:2) / Params.convert;
    Ta(:,3) = Ta(:,3) / (params.d_list(2)-params.d_list(1));

    % Add interpolation here
    Ta = round(Ta);

    mat_sz = max(Ta) - min(Ta) + 1;
    Tip_atoms2 = Ta - min(Ta) + 1;
    Tip_mat = zeros(mat_sz);
    for i = 1:size(Ta,1)
        Tip_mat(Tip_atoms2(i,1),Tip_atoms2(i,2),Tip_atoms2(i,3)) = exp(-params.Decay * d(i));
    end

    Force_mat2 = convn(Force_mat,flip(Tip_mat,3),'full');
    Force_mat2 = Force_mat2(mat_sz(1):end-mat_sz(1)+1,...
                            mat_sz(2):end-mat_sz(2)+1,...
                            mat_sz(3):end);
    z_mat = FindZeros(params.d_list, Force_mat2, params.init);
    I = z_mat';
    
    I = ZmapAdjust(I, params.adatom_offset, params.vacancy_offset,...
                      params.adatom_sigma, params.vacancy_sigma);
    
    %% Add noise
    sigma_1 = params.Sigma1;
    offset_xy_prob = params.drift_prob;
    offset_xy = params.drift_dist; % px
    offset_xy(1) = round(offset_xy(1));
    gauss1d_win = round(params.Gauss_win);
    gauss1d_sigma = params.Gauss_sig;
    gauss2d = params.Gauss2d;
    sigma_2 = params.Sigma2; % *sigma_1
    autoI = params.Autoint;
    if autoI == 0; autoI = 1e-5; end

    %% Parameter conversion, initiation, etc.
    sigma_2 = sigma_2 * sigma_1;
    padding = gauss1d_win;

    % Calibrate peak intensity
    Im = imgaussfilt(I,2);
    mM = mode(round(Im(:)*100))/100; % round to 0.01

    % Gaussian noise
    I = normrnd(I,sigma_1);

    % 1D scan texture
    w = gausswin(gauss1d_win , gauss1d_sigma);
    w = w/sum(w);

    for i = 1:size(I,1) % each row
        % random offset
        rand_x = 0; rand_y = 0;
        if rand() < offset_xy_prob(1)
            rand_x = randi([-1,1]*offset_xy(1));
        end
        if rand() < offset_xy_prob(2)
            rand_y = normrnd(0,offset_xy(2));
            % clip to [-1,1]
            if rand_y>1; rand_y = 1; end
            if rand_y<-1; rand_y=-1; end
        end

        % apply offset
        if rand_y~=0
            if (i>1 && rand_y<0) || i==size(I,1)
                temp2 = padarray(I(i-1,:)', padding, 'symmetric');
            else
                temp2 = padarray(I(i+1,:)', padding, 'symmetric');
            end            
        end

        % blur along scanning direction
        temp = padarray(I(i,:)', padding, 'symmetric');
        if rand_y~=0
            temp = temp * (1-abs(rand_y)) +...
                   temp2 * abs(rand_y);
        end
        temp = filter(w, 1, temp);

        sz_temp = length(temp);
        sz_pad_temp = floor([1.5,-0.5] * padding);
        I(i,:) = temp( 1 +sz_pad_temp(1)+ rand_x  : sz_temp +sz_pad_temp(2)+ rand_x  );
    end

    % Gaussian blur
    if gauss2d>0
        I = imgaussfilt(I,gauss2d);
    end
    % Gaussian noise
    I = normrnd(I,sigma_2);

%     I(I < (mM-autoI)) = mM-autoI;
%     I(I > (mM+autoI)) = mM+autoI;
    cm = mM - autoI;
    cM = mM + autoI;
    
    int_temp = (params.IntShift-0.5)*(cM-cm);
    cm = cm + int_temp;
    cM = cM + int_temp;
    
    %% Random contrast/brightness enhance/weaken patches
    numPatches = params.numPatches; % Adjust as needed
    blendingRadius = params.blendingRadius; % Blending region around patches
    
    originalImage = I;
    for i = 1:numPatches
        % Randomly select the center and radii for the oval patch
        [rows, cols] = size(originalImage);
        centerX = params.patchParams(i,1);
        centerY = params.patchParams(i,2);
        majorRadius = params.patchParams(i,3);
        minorRadius = params.patchParams(i,4);

        % Randomly select a contrast enhancement factor
        contrastFactor = params.patchParams(i,5);
        brightnessFactor = params.patchParams(i,6);

        % Create a binary mask for the oval patch
        [X, Y] = meshgrid(1:cols, 1:rows);
        ovalMask = ((X - centerX) / majorRadius).^2 + ((Y - centerY) / minorRadius).^2 <= 1;

        % Extract the patch from the original image
        patch = originalImage .* ovalMask;

        % Apply the contrast enhancement to the patch
        mI = mean(patch(patch>0));
        sI = std(patch(patch>0));
        enhancedPatch = (originalImage-mI) * contrastFactor + mI + brightnessFactor*sI;

        % Create an alpha mask for blending
        alphaMask = double(ovalMask);
        alphaMask = imgaussfilt(alphaMask, blendingRadius);
        alphaMask = alphaMask / max(alphaMask(:));
        
        % Apply alpha blending to smoothly transition the patch into the original image
        I = I .* (1 - alphaMask) + enhancedPatch .* alphaMask;
    end
end

function I = GTGeneration(offset)
    global atom_pos Params
    
    sel0 = atom_pos(:,1)==1 & abs(atom_pos(:,4))<Params.unmoved_Al;
    sel1 = atom_pos(:,1)==1 & atom_pos(:,4)>Params.unmoved_Al;
    
    
    % Scatter Al atom position, checking.
    if Params.scatter_check
        sel = sel0;
        figure()
%         scatter3(atom_pos(sel,2), atom_pos(sel,3), atom_pos(sel,4),'.')
%         scatter(atom_pos(sel,2), atom_pos(sel,3),'.')
        scatter(atom_pos(sel,2)/Params.convert + offset(1), ...
                atom_pos(sel,3)/Params.convert + offset(2),'.')
        axis equal
        axis([0 Params.generation_size 0 Params.generation_size])
    end
    
    % Create image
    I = zeros(Params.generation_size,Params.generation_size,3);
    [x, y] = meshgrid(1:Params.generation_size, 1:Params.generation_size);
    
    % Define Gaussian features
    numFeatures = sum(sel0);
    sel = find(sel0);
    for i = 1:numFeatures
        % Calculate the Gaussian intensity
        intensity = Params.Al_intensity * ...
                    exp(-((x - atom_pos(sel(i),2)/Params.convert - offset(1)).^2 +...
                          (y - atom_pos(sel(i),3)/Params.convert - offset(2)).^2) /...
                    (2 * Params.Al_sigma^2));
        % Add the Gaussian feature to the image
        I(:,:,1) = I(:,:,1) + intensity;
    end
    numFeatures = sum(sel1);
    sel = find(sel1);
    for i = 1:numFeatures
        % Calculate the Gaussian intensity
        intensity = Params.Al_intensity * ...
                    exp(-((x - atom_pos(sel(i),2)/Params.convert - offset(1)).^2 +...
                          (y - atom_pos(sel(i),3)/Params.convert - offset(2)).^2) /...
                    (2 * Params.Al_sigma^2));
        % Add the Gaussian feature to the image
        I(:,:,2) = I(:,:,2) + intensity;
    end
end

function Tip_atoms = GenerateTipAtoms(params)

    Tip_atoms = generate_Si(params.Tip_radius, params.Tip_height, 8, params.Facet111);

    atomic_positions = Tip_atoms;

    d = sqrt(sum(atomic_positions(:,1:2).^2,2));
    atomic_positions(:,3) = atomic_positions(:,3) + d.^2 * params.SphericalTip;
    atomic_positions(atomic_positions(:,3) > params.Tip_height, :) = [];
    params.Tip_atoms = atomic_positions;
end

function n = uniformRandomNumber(m,M)
    n = m + rand() * (M-m);
end

function n = normalRandomNumber(mu,sig)
    n = mu + randn * sig;
end



