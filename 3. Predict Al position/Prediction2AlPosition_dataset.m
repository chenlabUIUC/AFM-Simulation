%% Prediction2AlPosition_dataset.m
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

%% Input
% Real prediction
% h5file = 'prediction_AFM data 1.h5';
h5file = 'prediction_validation_20240312.h5';
data = h5read(h5file, '/data');
data = permute(data, [3,2,1,4]);
data = data(:,:,[2,3,1],:);
shape = size(data);

% Prepare output folder
output_dir = ['output/'];
CheckCreateDir(output_dir)
output_dir1 = ['output/1.Raw+Prediction/'];
CheckCreateDir(output_dir1)
output_dir2 = ['output/2.Upsample/'];
CheckCreateDir(output_dir2)
output_dir3 = ['output/3.Local maxima/'];
CheckCreateDir(output_dir3)

%% Preparation
% Generate patch series, Guassian peak with offset
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

figure(5);
set(gcf,'Position',[100,100,900,900]);
set(gca,'Position',[0 0 1 1]);

% for img_ind = 1:shape(4)
for img_ind = 15 % One image as example
    %% Main
    I = data(:,:,1,img_ind);
    I_upsample_base = imageUpsample(I, Patches, shape);
    isLocalMaxima_base = regionalMaxima(I_upsample_base);
    I = data(:,:,2,img_ind);
    I_upsample_adatom = imageUpsample(I, Patches, shape);
    isLocalMaxima_adatom = regionalMaxima(I_upsample_adatom);

    %% Output
    % Raw + Prediction
    imwrite(data(:,:,:,img_ind),[output_dir1,num2str(img_ind),'.tif']);
    
    % Up sampling
    indexedImage = ind2rgb(gray2ind(I_upsample_base, 256), jet(256));
    rgbImage = uint8(indexedImage * 255);
    imwrite(rgbImage, [output_dir2,num2str(img_ind),'-base.tif']);
    indexedImage = ind2rgb(gray2ind(I_upsample_adatom, 256), jet(256));
    rgbImage = uint8(indexedImage * 255);
    imwrite(rgbImage, [output_dir2,num2str(img_ind),'-adatom.tif']);

    % Local maxima
    rec = [];
    hold off; imshow(I_upsample_base); hold on;
    [rows, cols] = find(isLocalMaxima_base & I_upsample_base > Params.Intensity_thresh);
    scatter(cols,rows,'yx')
    saveas(gcf,[output_dir3,num2str(img_ind),'-base.tif'])
    rec = cat(1,rec,cat(2,cols, rows, ones(size(cols))));
    
    hold off; imshow(I_upsample_adatom); hold on;
    [rows, cols] = find(isLocalMaxima_adatom & I_upsample_adatom > Params.Intensity_thresh);
    scatter(cols,rows,'yx')
    saveas(gcf,[output_dir3,num2str(img_ind),'-adatom.tif'])
    rec = cat(1,rec,cat(2,cols, rows, ones(size(cols))*2));
    
    outputFileName = [output_dir3,num2str(img_ind),'.txt'];
    fileID = fopen(outputFileName, 'w');
    fprintf(fileID, '%d %d %d\n', rec');
    fclose(fileID);
end

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

