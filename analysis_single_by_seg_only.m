% By Ashkan Pakzad 25th July 2019

% attempt to identify bifurcations from global segmentation and skel2graph
% alone.

%% load tapering data & Segmentation
load('P_N1_raw_image_CT_nifti_A_22195509.mat')
segname = 'N1_raw_image_lumen_seg.nii';
S = logical(niftiread(segname));

%% convert seg to graph
skel = bwskel(S);

[A,node,link] = Skel2Graph3D(skel,0); % dependent function
G = graph(A);

%% identify carina and terminal point node.
% last node corresponds to top of trachea.
carina_idx = neighbors(G,length(A)); % carina node
all_term_idx = find([node.ep]);
terminal_vox = [node(all_term_idx).idx]; % vox indicies of terminal points

terminal_voxidx = find(terminal_vox == sturct_to_save.airway_number); % the idx for desired terminal
terminal_idx = all_term_idx(terminal_voxidx);

%% identify path of airway
[P, d] = shortestpath(G,carina_idx,terminal_idx);
disp(d)
disp(sturct_to_save.tapering_raw_image.arclegth(end))
bifurcations = length(P) - 2; % first and last are terminal. 
generation = bifurcations + 1;

dist = zeros(length(P)-1, 1);
for i = 1:length(P)-1
    if i == 1
        [~, dist(i)] = shortestpath(G,P(i),P(i+1));
    else
        [~, dd] = shortestpath(G,P(i),P(i+1));
        dist(i) = dd+dist(i-1);
    end
end

%dist = dist + 10;

%% visualise bifucation points against area from tapering
figure
plot(sturct_to_save.tapering_raw_image.arclegth,...
    sturct_to_save.tapering_raw_image.area_results.phyiscal_area)
hold on
plot(dist, 20, 'r.', 'MarkerSize', 14)


% OBSERVED VISUALLY AT: 67.6, 79.8, 84.45, 105.6, 115.8, 128.4, 132.7,
% 143,159, 181(?) - 9/10 bifurcations

% SHORTEST PATH:                      [57;67;82;91;101;110;116;125;142;180]
% SHORTEST PATH + 10 (CARINA start):  [67;77;92;101;111;120;126;135;152;190]

% figure
% yy = smooth(lumen_major);
% plot(arc_length, lumen_major);
% hold on
% plot(arc_length, yy);
%% split taper data into segments

%% Find path voxels
voxel_in_path = [];
for i = 1:length(P)-1
    current_link = find([link.n1]== P(i+1) & [link.n2]== P(i) |...
        [link.n1]== P(i) & [link.n2]== P(i+1));
    voxel_in_path = cat(1,voxel_in_path, [link(current_link).point]');
end

[J,I,K] = ind2sub(size(skel), voxel_in_path); % flipped from skeleton for some reason

%% Visualise path from graph on skeleton
%path3D = zeros(size(skel));
%path3D(voxel_in_path) = 1;

[strstopJ, strstopI, strstopK] = ind2sub(size(skel),[node(carina_idx).idx, sturct_to_save.airway_number]);

figure
isosurface(skel)
hold on
plot3(I,J,K,'r.', 'MarkerSize',15)
plot3(strstopI(1),strstopJ(1),strstopK(1),'mo','MarkerSize',20)
plot3(strstopI(end),strstopJ(end),strstopK(end),'mo','MarkerSize',20)
view(90,0)



%% extract
arc_length = sturct_to_save.tapering_raw_image.arclegth;
area_info = sturct_to_save.tapering_raw_image.area_results;
lumen_area = [area_info.phyiscal_area];
wall_area = [area_info.phyiscal_area_wall];

% ellipse, stored as (Cx, Cy, Rx, Ry, theta_radians)
lumen_ellipse = {area_info.elliptical_info{:,1}}';
wall_ellipse = {area_info.elliptical_info_wall{:,1}}';
lumen_radius = zeros(2, numel(arc_length));
wall_radius = zeros(2, numel(arc_length));
for i = 1:numel(arc_length)
    lumen_radius(1,i) = lumen_ellipse{i, 1}.elipllical_info(3);
    lumen_radius(2,i) = lumen_ellipse{i, 1}.elipllical_info(4);
    wall_radius(1,i) = wall_ellipse{i, 1}.elipllical_info(3);
    wall_radius(2,i) = wall_ellipse{i, 1}.elipllical_info(4);
end
lumen_major = max(lumen_radius);
wall_major = max(wall_radius);

