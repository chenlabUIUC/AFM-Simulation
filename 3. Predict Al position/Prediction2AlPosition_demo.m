%% Prediction2AlPosition_demo.m
% Process prediction to get Al positions, using one frame from the
% training/validation set.
%
% INPUT:
%   Input files: 'prediction_xxx.h5'
%   Parameters of how ground truth is generated:
%       'Params.Al_intensity', 'Params.Al_sigma', 'Params.convert' as
%       described in dataset generation part.
%   Parameters related to upsampling:
%       'Params.Upsampling': how many times of resolution upsampling.
%       'Parmas.PatchSize': the patch size for upsampling.
%   Interpreting the prediction result:
%       'Params.Intensity_thresh': larger than which is identified as
%       atoms.
%   For demonstrating the performance, more GT and random parameters are
%   needed: 'GT_validation set.mat', 'Params_training set+validation
%   set.mat'.
%
% OUTPUT:
%   All output images are in the created 'output/' folder.
%
% HISTORY:  Written by Chang (Chris) Qian Last modified by Chang (Chris)
% Qian on 04/04/2024



clear; close all; clc;

global Params
Params.Al_intensity = 1;
Params.Al_sigma = 0.5;
Params.Upsampling = 8;
Params.PatchSize = 9;
Params.Intensity_thresh = 0.15;
Params.convert = 200/256; % A/px

% Check on validation set
h5file = 'prediction_validation_20240312.h5';
data = h5read(h5file, '/data');
data = permute(data, [3,2,1,4]);
data = data(:,:,[2,3,1],:);
shape = size(data);

% Load GT
load('GT_validation set.mat');
load('Params_training set+validation set.mat');

% Generate patch series
Patches = zeros(Params.PatchSize, Params.PatchSize, Params.Upsampling^2);
[x, y] = meshgrid(1:Params.PatchSize, 1:Params.PatchSize);
for x_i = 1:Params.Upsampling
    for y_i = 1:Params.Upsampling
        atom_pos = [0.5, 0.5]*Params.PatchSize - [0.5,0.5] + [x_i-1,y_i-1] / Params.Upsampling;
        ind = (x_i-1)*Params.Upsampling + y_i;
        intensity = Params.Al_intensity * ...
                    exp(-((x - atom_pos(1)).^2 +...
                         (y - atom_pos(2)).^2) /...
                    (2 * Params.Al_sigma^2));
        Patches(:,:,ind) = Patches(:,:,ind) + intensity;
    end
end

%% Demo
figure(); tiledlayout(Params.Upsampling,Params.Upsampling, 'TileSpacing','tight','Padding','tight')
set(gcf,'Position',[100,100,900,900]);
for i = 1:size(Patches,3)
    nexttile()
    row = floor((i-1)/Params.Upsampling);
    col = mod((i-1) , Params.Upsampling);
    ind = col*Params.Upsampling + row + 1;
    imshow(Patches(:,:,ind))
end


for img_ind = 15
    I = data(:,:,1,img_ind);
    I_upsample_base = imageUpsample(I, Patches, shape);
    isLocalMaxima_base = regionalMaxima(I_upsample_base);
    I = data(:,:,2,img_ind);
    I_upsample_adatom = imageUpsample(I, Patches, shape);
    isLocalMaxima_adatom = regionalMaxima(I_upsample_adatom);


    zoom = [100.0623  115.7866   51.4214   67.1457];

    figure(); tiledlayout(3,3,"TileSpacing","tight","Padding","compact");
    set(gcf,'Position',[100,100,900,900]);
    nexttile(); imshow(data(:,:,3,img_ind));
    nexttile(); 
        I = GT(img_ind).I; I(:,:,3) = 0;
        I = flipud(I);
        imshow(I);
        drawRectangle(zoom)
    nexttile();
        imshow(I);
        axis(zoom)
        hold on;

        mat_sz = getMatSz(params_rec(img_ind + 500));
        pos = GT(img_ind).pos;
        pos = pos/Params.convert;
        offset1 = flip(params_rec(img_ind + 500).offset_I);
        offset2 = [-floor(mat_sz(1)/2)+0.5, -floor(mat_sz(2)/2)+0.5];
        pos(:,[1:2]) = pos(:,[1:2]) - offset1 + offset2;
        pos(:,2) = size(I,1) - pos(:,2) +1;
        % pos(:,[1:2]) = pos(:,[1:2]) + [5,5] + [-1,0];
        

        select= pos(:,3)==0;
        scatter(pos(select,1),pos(select,2),'yx','LineWidth',2)
        select= pos(:,3)>0;
        scatter(pos(select,1),pos(select,2),'cx','LineWidth',2)
    nexttile(); axis off;
    nexttile();
        imshow(makergb(data(:,:,[1:2],img_ind),[1,2]))
        drawRectangle(zoom)
    nexttile();
        imshow(makergb(data(:,:,[1:2],img_ind),[1,2]))
        axis(zoom)
    nexttile(); axis off;
    nexttile();
        imshow(makergb(cat(3,I_upsample_base,I_upsample_adatom),[1,2]))
        drawRectangle(zoom.*Params.Upsampling - Params.PatchSize)
    nexttile();
        imshow(makergb(cat(3,I_upsample_base,I_upsample_adatom),[1,2]))
        axis(zoom.*Params.Upsampling - Params.PatchSize)
        hold on
        [rows, cols] = find(isLocalMaxima_base & I_upsample_base > Params.Intensity_thresh);
        scatter(cols,rows,'yx','LineWidth',2)
        [rows, cols] = find(isLocalMaxima_adatom & I_upsample_adatom > Params.Intensity_thresh);
        scatter(cols,rows,'cx','LineWidth',2)
end

%% Error estimation
error_base = [];
error_adatom = [];
for img_ind = 1:20
    % Prediction Al atom position
    I = data(:,:,1,img_ind);
    I_upsample_base = imageUpsample(I, Patches, shape);
    isLocalMaxima_base = regionalMaxima(I_upsample_base);
    I = data(:,:,2,img_ind);
    I_upsample_adatom = imageUpsample(I, Patches, shape);
    isLocalMaxima_adatom = regionalMaxima(I_upsample_adatom);

    % GT Al atom position
    mat_sz = getMatSz(params_rec(img_ind + 500));
    pos = GT(img_ind).pos;
    pos = pos/Params.convert;
    offset1 = flip(params_rec(img_ind + 500).offset_I);
    offset2 = [-floor(mat_sz(1)/2)+0.5, -floor(mat_sz(2)/2)+0.5];
    pos(:,[1:2]) = pos(:,[1:2]) - offset1 + offset2;
    pos(:,2) = size(I,1) - pos(:,2) +1;

    [pos_base_y, pos_base_x] = find(isLocalMaxima_base & I_upsample_base > Params.Intensity_thresh);
    [pos_adatom_y, pos_adatom_x] = find(isLocalMaxima_adatom & I_upsample_adatom > Params.Intensity_thresh);

    pos_test = [pos_base_x, pos_base_y];
        pos_test = pos_test + Params.PatchSize - [1,1];
        pos_test = pos_test ./ Params.Upsampling;
        for i = 1:size(pos_test,1)
            select = pos(:,3)==0;
            e = pos_test(i,:) - pos(select,1:2);
            ed = sqrt(sum(e.^2,2));
            [em,b] = min(ed);
            error_base = [error_base; e(b,:), em];
        end

    pos_test = [pos_adatom_x, pos_adatom_y];
        pos_test = pos_test + Params.PatchSize - [1,1];
        pos_test = pos_test ./ Params.Upsampling;
        for i = 1:size(pos_test,1)
            select = pos(:,3)>0;
            e = pos_test(i,:) - pos(select,1:2);
            ed = sqrt(sum(e.^2,2));
            [em,b] = min(ed);
            error_adatom = [error_adatom; e(b,:), em];
        end
end

figure(); tiledlayout(1,3,'TileSpacing','tight','Padding','compact');
set(gcf,'Position',[100,100 1200 400])
nexttile();
    scatterhistogram(error_base(:,1), error_base(:,2), ...
        'XLimits',[-2,2],'YLimits',[-2,2], ...
        'HistogramDisplayStyle','smooth','MarkerStyle','.','Color','r')
nexttile();
    scatterhistogram(error_adatom(:,1), error_adatom(:,2), ...
        'XLimits',[-5,5],'YLimits',[-5,5], ...
        'HistogramDisplayStyle','smooth','MarkerStyle','.','Color','g')
nexttile();
    hold on;
    histogram(error_base(:,3),'BinLimits',[0,3],'Normalization','pdf','LineStyle','none','FaceColor','r')
    histogram(error_adatom(:,3),'BinLimits',[0,7],'Normalization','pdf','LineStyle','none','FaceColor','g')

disp(['Error for base Al atoms： ', num2str(mean(error_base(error_base(:,3)<3,3))),' nm'])
disp(['Error for Al adatoms: ', num2str(mean(error_adatom(error_adatom(:,3)<7,3))), ' nm'])


%% Functions
function I_upsample = imageUpsample(I, Patches, shape)
    global Params
    I_upsample = zeros(shape(1:2) * Params.Upsampling);
    for x_i = 1:Params.Upsampling
        for y_i = 1:Params.Upsampling
            ind = (x_i-1)*Params.Upsampling + y_i;
            I_upsample((Params.Upsampling-y_i+1):Params.Upsampling:end,...
                       (Params.Upsampling-x_i+1):Params.Upsampling:end) = ...
                conv2(I,Patches(:,:,ind),'same');
        end
    end
end

function I = regionalMaxima(I_upsample)
    global Params

    originalImage = I_upsample;

    % Convert the image to grayscale if it's not already in grayscale
    if size(originalImage, 3) == 3
        originalImage = rgb2gray(originalImage);
    end

    % Define the half-width and half-height of the rectangular neighborhood
    halfWidth = Params.Upsampling;
    halfHeight = Params.Upsampling;

    % Initialize a binary image to store the local maxima
    isLocalMaxima = zeros(size(originalImage));

    % Get the size of the image
    [rows, cols] = size(originalImage);

    % Iterate through each pixel in the image
    for i = halfHeight+1 : rows-halfHeight
        for j = halfWidth+1 : cols-halfWidth
            if originalImage(i, j) < Params.Intensity_thresh
                continue
            end
            
            % Extract the rectangular neighborhood
            neighborhood = originalImage(i-halfHeight:i+halfHeight, j-halfWidth:j+halfWidth);

            % Check if the pixel value is the maximum in the neighborhood
            if originalImage(i, j) == max(neighborhood(:))
                isLocalMaxima(i, j) = 1;
            end
        end
    end
    
    I = isLocalMaxima;
end

function CheckCreateDir(dir)
    if ~exist(dir,'dir'); mkdir(dir); end
end

function I = makergb(I0,ax)
    I = zeros([size(I0,[1,2]),3]);
    if size(I0,3) == 1
        I(:,:,ax) = I0;
    else
        for i = 1:length(ax)
            I(:,:,ax(i)) = I0(:,:,i);
        end
    end
end

function drawRectangle(zoom)
    pos = [zoom([1,3]), zoom(2)-zoom(1), zoom(4)-zoom(3)];
    rectangle('Position', pos, 'EdgeColor','w','LineWidth',2);
end

function m = getMatSz(params)
rot_mat = @(x) [cosd(x), -sind(x), 0;
                sind(x), cosd(x),  0;
                0,       0,        1;];
convert = 200/256;
    Ta = params.Tip_atoms * rot_mat(params.Rot);
    Ta(:,1:2) = Ta(:,1:2) / convert;
    Ta(:,3) = Ta(:,3) / (params.d_list(2)-params.d_list(1));
    Ta = round(Ta);
    m = max(Ta) - min(Ta) + 1;
end


