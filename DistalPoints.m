% By Ashkan Pakzad 4th July 2019 (ashkan.pakzad.13 (at) ucl.ac.uk)

% Marks the distal terminals of airways and save as nii for
% provided airway segmentation. 
% Assumes the trachea is fully segmented (beyond the lungs) and at the top 
% of the segmentation image.

% DEPENDANCY:
%   Skel2Graph3D: https://github.com/phi-max/skel2graph3d-matlab

% INPUT 
%   segname: segmentation file name
%   distalname: desired name for new distal array saved as .nii
%   branch_threshold: for input to skel2graph3D, set the minimum branch
%   length that should be considered. (used to prevent false branches)
%   n: number of terminal points to consider. 
%       (set low for testing, otherwise set 0 for all airways.)
% OUTPUT
%   binary nifti image with terminals labelled as 1.

% EXAMPLE USE:
%   create a distal.nii image based on airway segmentation called
%   airway_seg.nii. Consider every branch and only return the first three
%   distal points.
% DistalPoints('airway_seg.nii','distal',0,3)

function DistalPoints(segname,distalname,branch_threshold,n)

%% Load segmentation
S = logical(niftiread(segname));

%% Skeletonise lumen segmentation & graph it
skel = bwskel(S);

[~,node,~] = Skel2Graph3D(skel,branch_threshold); % dependent function

%% identify terminal points
terminal_idx = find([node.ep]);

% remove trachea
[~, trachea_idx] = min([node(terminal_idx).comz]);
terminal_idx(terminal_idx == trachea_idx) = 0;
terminal_idx = terminal_idx(1:end-1);

if n ~= 0
    terminal_idx = terminal_idx(1:n);
else
end

%% Convert to voxel list and create logical array
try
    terminal_vox = [node(terminal_idx).idx];
catch
    warning('Terminal nodes appear to be assigned to multiple voxels.')
end

[i,j,k] = ind2sub(size(S),terminal_vox);

flip_terms = sub2ind(size(S),i,j,k);

distal = zeros(size(S));
distal(flip_terms) = 1;

%% Save distal logical array
niftiwrite(double(distal), distalname)
