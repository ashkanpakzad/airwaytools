% By Ashkan Pakzad 25th July 2019

% attempt to identify bifurcations from tapering segmentation.

%% load tapering data & Segmentation
load('P_N1_raw_image_CT_nifti_A_22195509.mat')
segname = 'N1_raw_image_lumen_seg.nii';
S = logical(niftiread(segname)); % NOTE upper and lower s, S
s_raw = sturct_to_save.tapering_seg_image.resample_image;

arc_length = sturct_to_save.tapering_raw_image.arclegth;
area_info = sturct_to_save.tapering_raw_image.area_results;
lumen_area = [area_info.phyiscal_area];
wall_area = [area_info.phyiscal_area_wall];

%% remove bronchi at carina
removal = 50;
s_raw = s_raw(:,:,removal:end);

%% Preprocessing of indiv. airway segmentation
% force binary image
sraw_binary = s_raw > 0.7;

% smoothen binary image
N = 3;
kernel = ones(N, N, N) / N^3;
blurryImage = convn(double(sraw_binary), kernel, 'same');
s_binary = blurryImage > 0.5;

% fill holes
s_fill = imfill(s_binary, 'holes');

% remove smaller structures.
s = zeros(size(s_fill));
CC = bwconncomp(s_fill);
numPixels = cellfun(@numel,CC.PixelIdxList); % apply func to each cell
[~,idx] = max(numPixels);
s(CC.PixelIdxList{idx}) = 1;
s = logical(s);

%% skeletonise and convert single airway to graph
s_skel = bwskel(s);

[sA,snode,slink] = Skel2Graph3D(s_skel,0); % dependent function
sG = graph(sA);

%% visualise airway skeleton
SE = strel('disk', 1);
s_skel_display = imdilate(s_skel,SE);

figure
plotGraphSkel(s_skel_display, snode) 

%% identify start and terminal node
s_all_term_idx = find([snode.ep]);
sterminal_vox = [snode(s_all_term_idx).idx]; % vox indicies of terminal points
[sJ,sI,sK] = ind2sub(size(s_skel),sterminal_vox);

% first z node = start
[~,start_zidx] = min(sK');
start_idx = s_all_term_idx(start_zidx);

% last z node = airway end
[~,stop_zidx] = max(sK');
stop_idx = s_all_term_idx(stop_zidx);

sP = shortestpath(sG,start_idx,stop_idx); % path for airway

%% remove false bifurcations by analysing skeleton
% Func. to find plausible bifurcations
[true_bifurcations, ~] = Bifurcation_Plausibility(sP, snode, slink);
% extract plausible bifurcations
sP_bi = sP(2:end-1); % removes start and end nodes
sP_plaus = sP_bi(true_bifurcations);
% re-add start and end nodes
sP_plaus = [sP(1), sP_plaus, sP(end)];

% counter for false bifurcations
false_counter = sum(~true_bifurcations);

%% identify arc_length distances for path bifurcation nodes
path_zind = round([snode(sP_plaus).comz]);
%[pathJ,pathI,pathK] = ind2sub(size(s_skel),path_ind);
path_zind = path_zind(2:end-1); % remove start and end
bifurcation_idx = path_zind+removal;
bifurcation_points = arc_length(bifurcation_idx);

%% identify bifurcations with local maxima peaks raw
[raw_pks,raw_locs] = findpeaks(lumen_area);
new_bifurcation_idx = zeros(size(bifurcation_idx));
for i = 1:length(bifurcation_idx)
    current_idx = bifurcation_idx(i);
    % identify nearest peak
    [min_val, min_idx] = min(abs(raw_locs - bifurcation_idx(i)));
    try
        % consider if nearer peaks are more likely to rep. bifurcation
        [~, max_idx] = max(raw_pks(min_idx-3:min_idx+3));
        if max_idx > 4
            idx = min_idx - 4 + max_idx;
        elseif max_idx < 4
            idx = min_idx - (4 - max_idx);
        else
            idx = min_idx;
        end
    catch
         idx = min_idx;
    end
    raw_pks(idx) = 0; % prevents same peak being chosen twice.
    new_bifurcation_idx(i) = raw_locs(idx);
    % consider something here incase finding peak fails
end

%% visualise all possible peaks and original idx
figure
plot(lumen_area)
hold on
plot(bifurcation_idx, lumen_area(bifurcation_idx), 'c.', 'MarkerSize', 14)
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
findpeaks(lumen_area)

%% visualise new peak bifurcations and original idx
figure
plot(arc_length,lumen_area)
hold on
plot(bifurcation_points, lumen_area(bifurcation_idx), 'r.', 'MarkerSize', 14)
plot(arc_length(new_bifurcation_idx),lumen_area(new_bifurcation_idx), '^')
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Raw Area','new new Bifurcation Points')

%% Remove Bifurcation Points
% most computationally intensive.
% add carina to be excluded.
new_bifurcation_idx = [1 new_bifurcation_idx];
%indices for exclusion between lower and upper point at each bifurcation.
lower_exclude = zeros(size(new_bifurcation_idx));
upper_exclude = zeros(size(new_bifurcation_idx));

for i = 1:length(new_bifurcation_idx)
    
    pck_idx = new_bifurcation_idx(i);
    try
        if i == 1
            data_hwidth = 100;% initial peak width for carina
        else
            data_hwidth = 50;% initial peak width
        end
        % initiate loop variables
        gof_n = 0;
        gof_np1 = 0.0001;
        int = 1; % inteval width reduction
        lower_width = pck_idx-data_hwidth;
        upper_width = pck_idx+data_hwidth;
        % set precautions incase 1 < data width > last idx.
        if lower_width < 1
            lower_width = 1;
        else
        end
        if upper_width > length(lumen_area)
            upper_width = length(lumen_area);
        else
        end
        x_data = mat2gray(lumen_area(lower_width:upper_width));
        [~, x_data_pck_idx] = max(x_data);
        
         % reduce data window if bifurcation is not the peak.
        if i ~= 1
            while x_data_pck_idx ~= data_hwidth+1
                [~, x_data_pck_idx] = max(x_data);
                data_hwidth = data_hwidth - int;
            end
        else
        end

        while gof_n < gof_np1 || gof_n < 0.2
            try % save prev. result for display.
                n_fitresult = fitresult; 
                n_data_hwidth = data_hwidth;
                n_x_data = x_data;
                n_vals = vals;
                catch
            end
            
            % normalise data window
            if lower_width < 1
                lower_width = 1;
            else
            end
            if upper_width > length(lumen_area)
                upper_width = length(lumen_area);
            else
            end
            x_data = mat2gray(lumen_area(lower_width:upper_width));
            
            % equally spaced values
            vals = [1:length(x_data)]';

            % fit gaussian with fixed mean at centre
            options = fitoptions('gauss1');
            if i == 1 % mean is 1 for carina
                options.Lower = [-inf 1 0];
                options.Upper = [inf 1 inf];
            else
                options.Lower = [-inf data_hwidth+1 0];
                options.Upper = [inf data_hwidth+1 inf];
            end
            [fitresult, gof] = fit(vals, x_data, 'gauss1', options);
            
            % assign values for next loop
            gof_n = gof_np1; 
            gof_np1 = gof.rsquare;
            % reduce peak width 
            data_hwidth = data_hwidth - int;
            lower_width = pck_idx-data_hwidth;
            upper_width = pck_idx+data_hwidth;
        end
    catch
        n_data_hwidth = 3; % if fail to fit, remove 7 points around bifurcation
    end
    if i == 1
        lower_exclude(i) = 1;
    else
        lower_exclude(i) = pck_idx-n_data_hwidth;
    end
    upper_exclude(i) = pck_idx+n_data_hwidth;
    figure
    plot(n_vals,n_x_data,'.')
    hold on
    plot(n_fitresult)
    disp(['Goodness of fit, rsquare = ', num2str(gof.rsquare)])
end

% figure
% plot(n_vals,n_x_data,'.')
% hold on
% plot(n_fitresult)
% disp(['Goodness of fit, rsquare = ', num2str(gof.rsquare)])

%% extract all points not excluded
logic_include = logical(ones(size(lumen_area)));
for i = 1:length(lower_exclude)
    logic_include(lower_exclude(i):upper_exclude(i)) = 0;
end
bi_removed_arcl = arc_length(logic_include);
bi_removed_area = lumen_area(logic_include);
bi_removed_warea = wall_area(logic_include);

%% visualise excluded and included points
figure
subplot(2,1,1)
plot(arc_length, lumen_area, '.')
hold on
plot(bi_removed_arcl, bi_removed_area, '.')
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Points','Non-Bifurcation points')
title('lumen area')

subplot(2,1,2)
plot(arc_length, wall_area, '.')
hold on
plot(bi_removed_arcl, bi_removed_warea, '.')
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Points','Non-Bifurcation points')
title('WALL area')


% %% Create global airway skeleton and convert to graph
% S_skel = bwskel(S);
% 
% [SA,Snode,Slink] = Skel2Graph3D(S_skel,0); % dependent function
% SG = graph(SA);
% 
% %% identify carina and terminal point node.
% % last node corresponds to top of trachea.
% carina_idx = neighbors(SG,length(SA)); % carina node
% all_term_idx = find([Snode.ep]);
% terminal_vox = [Snode(all_term_idx).idx]; % vox indicies of terminal points
% 
% terminal_voxidx = find(terminal_vox == sturct_to_save.airway_number); % the idx for desired terminal
% terminal_idx = all_term_idx(terminal_voxidx);
% 
% %% identify path of airway
% [P, d] = shortestpath(SG,carina_idx,terminal_idx);
% disp(d)
% disp(sturct_to_save.tapering_raw_image.arclegth(end))
% bifurcations = length(P) - 2; % first and last are terminal. 
% generation = bifurcations + 1;
% 
% dist = zeros(length(P)-1, 1);
% for i = 1:length(P)-1
%     if i == 1
%         [~, dist(i)] = shortestpath(SG,P(i),P(i+1));
%     else
%         [~, dd] = shortestpath(SG,P(i),P(i+1));
%         dist(i) = dd+dist(i-1);
%     end
% end
% 
% %dist = dist + 10;
% 
% %% split taper data into segments
% 
% %% Find path voxels
% voxel_in_path = [];
% for i = 1:length(P)-1
%     current_link = find([Slink.n1]== P(i+1) & [Slink.n2]== P(i) |...
%         [Slink.n1]== P(i) & [Slink.n2]== P(i+1));
%     voxel_in_path = cat(1,voxel_in_path, [Slink(current_link).point]');
% end
% 
% [J,I,K] = ind2sub(size(S_skel), voxel_in_path); % flipped from skeleton for some reason
% 
% %% Visualise path from graph on skeleton
% %path3D = zeros(size(skel));
% %path3D(voxel_in_path) = 1;
% 
% [strstopJ, strstopI, strstopK] = ind2sub(size(S_skel),[Snode(carina_idx).idx, sturct_to_save.airway_number]);
% 
% figure
% isosurface(S_skel)
% hold on
% plot3(I,J,K,'r.', 'MarkerSize',15)
% plot3(strstopI(1),strstopJ(1),strstopK(1),'mo','MarkerSize',20)
% plot3(strstopI(end),strstopJ(end),strstopK(end),'mo','MarkerSize',20)
% view(90,0)
% 
% %% extract
% arc_length = sturct_to_save.tapering_raw_image.arclegth;
% area_info = sturct_to_save.tapering_raw_image.area_results;
% lumen_area = [area_info.phyiscal_area];
% wall_area = [area_info.phyiscal_area_wall];
% 
% % ellipse, stored as (Cx, Cy, Rx, Ry, theta_radians)
% lumen_ellipse = {area_info.elliptical_info{:,1}}';
% wall_ellipse = {area_info.elliptical_info_wall{:,1}}';
% lumen_radius = zeros(2, numel(arc_length));
% wall_radius = zeros(2, numel(arc_length));
% for i = 1:numel(arc_length)
%     lumen_radius(1,i) = lumen_ellipse{i, 1}.elipllical_info(3);
%     lumen_radius(2,i) = lumen_ellipse{i, 1}.elipllical_info(4);
%     wall_radius(1,i) = wall_ellipse{i, 1}.elipllical_info(3);
%     wall_radius(2,i) = wall_ellipse{i, 1}.elipllical_info(4);
% end
% lumen_major = max(lumen_radius);
% wall_major = max(wall_radius);
