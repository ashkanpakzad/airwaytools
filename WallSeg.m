% By Ashkan Pakzad on 9th July 2019

% Script to segment airway walls given an airway lumen input

%% Params
imagename='N1.nii.gz';
segname='N1_lumen';
r = 3; % = max wall thickness/2 in voxel.
% airway wall H.U. range.
LB = -800; UB = -100; 
% upper airway lumen H.U. range.
AB = -800;

%% load data & segmentations
V = double(niftiread(imagename));
L = logical(niftiread(segname));

%% Dialate segmentation
SE = strel('disk', r);
D = imdilate(L,SE);

%% Keep if part of airway.
T = V > LB & V < UB & D == 1;

%% combine lumen and threshold seg, fill holes & remove smaller structures
P = T == 1 | L == 1;
% 3D holes
F3 = imfill(P, 'holes');
% 2D holes
F2 = zeros(size(F3));
for i = 1:size(L,1)
    F2(i,:,:)= imfill(F3(i,:,:), 'holes');
end
for i = 1:size(L,2)
    F2(:,i,:)= imfill(F3(:,i,:), 'holes');
end
for i = 1:size(L,3)
    F2(:,:,i)= imfill(F3(:,:,i), 'holes');
end

% remove smaller structures
C = zeros(size(F2));
CC = bwconncomp(F2);
numPixels = cellfun(@numel,CC.PixelIdxList); % apply func to each cell
[~,idx] = max(numPixels);
C(CC.PixelIdxList{idx}) = 1;
C = logical(C);

%% reconsider lumen voxels
NL = V <= AB & C == 1 | L ==1;

%% extract wall
W = C - NL;

%% Visualise output
% figure
% isosurface(NL);
% hold on
% title('lumen')
% axis tight
% 
% figure
% isosurface(F); % use filled seg for visualisation
% hold on
% title('wall')
% axis tight

%% Save segmentation
niftiwrite(double(W), 'wall_v4')
niftiwrite(double(NL), 'lumen_v4')