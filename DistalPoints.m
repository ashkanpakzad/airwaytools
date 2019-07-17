% By Ashkan Pakzad 4th July 2019 (ashkan.pakzad.13@ucl.ac.uk)

% function to mark the distal terminals of airways and save as nii for
% provided segmentation name in .nii format.

% I - segmentation file name
%   - desired name for new distal array saved as .nii
%   - number of terminal points to use. 0 for all.

% example use:
% DistalPoints('N1_airway.nii.gz','test_distal',3)

function DistalPoints(segname,distalname,n)

%% Load data
S = logical(niftiread(segname));

%% Skeletonise lumen segmentation & graph it
skel = bwskel(S);

[A,node,link] = Skel2Graph3D(skel,0); % dependent function
G = graph(A);

%% Find carina
% last node APPARENTLY corresponds to top of trachea.
%C_Gidx = neighbors(G,length(A)); % carina node

%% identify terminal points
terminal_idx = find([node.ep]);

% remove trachea
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

%% visualise Graph plot

% figure
% G = graph(A);
% h = plot(G);
% highlight(h, C_Gidx); % highlight carina
% highlight(h, terminal_idx); % highlight carina

%% visualise terminal points
% [I,J,K]=ind2sub(size(distal),terminal_vox);
% 
% figure
% isosurface(S)
% hold on
% plot3(J,I,K,'.r','MarkerSize',16) % nb swap of axes for visualisation


