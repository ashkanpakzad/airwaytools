% By Ashkan Pakzad 4th July 2019 (ashkan.pakzad.13@ucl.ac.uk)

% Script to mark the distal terminals of airways.

%% Params
imagename='N1.nii.gz';
segname='N1_lumen.nii.gz';

%% Load data
V = double(niftiread(imagename));
I = niftiinfo(imagename);
S = logical(niftiread(segname));

%% Skeletonise lumen segmentation & graph it
skel = bwskel(S);

[A,node,link] = Skel2Graph3D(skel,0); % dependent function
G = graph(A);

%% Find carina
% last node APPARENTLY corresponds to top of trachea.
C_Gidx = neighbors(G,length(A)); % carina node


%% identify terminal points
terminal_idx = find([node.ep])';

% remove trachea
terminal_idx = terminal_idx(1:end-1);

%% Convert to voxel list and create logical array
terminal_vox = [node(terminal_idx).idx];

[I,J,K] = ind2sub(size(S), terminal_vox);

ind_idx = sub2ind(size(S),J,I,K);

distal = zeros(size(S));
distal(ind_idx) = 1; % indicies of I and J must be swapped to for CT space.

%% Save distal logical array
niftiwrite(double(distal), 'N1_distal')

%% tools to Graph plot

figure
G = graph(A);
h = plot(G);
highlight(h, C_Gidx); % highlight carina
highlight(h, terminal_idx); % highlight carina



