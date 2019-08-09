% Script to save N1 bi points.
% By Ashkan Pakzad on 30th July 2019

%% input var
data_file = fullfile('..','..','results',)
prefix = 'P_N1_raw_image_CT_nifti_A_';
% remove bronchi at carina
removal = 50;
search_area = 4; % for local maxima
airway_number = [N1_cleaned_data(:).airway_number];

%N1_cleaned_data(1).bi_loc = [];

%% run loop to determine bifurcation data points

for j = 52:133
    disp(['airway ', num2str(j) ' of ' , num2str(length(airway_number))])
    try % incase airway missing
    %% load
    filename = [prefix,num2str(airway_number(j)),'.mat'];
    load(filename);
    s_raw = sturct_to_save.tapering_seg_image.resample_image;
    
    % extract info from airway struct.
    arc_length = sturct_to_save.tapering_raw_image.arclegth;
    area_info = sturct_to_save.tapering_raw_image.area_results;
    lumen_area = [area_info.phyiscal_area];
    wall_area = [area_info.phyiscal_area_wall];
    
    %% Preprocessing of indiv. airway segmentation
    s_raw = s_raw(:,:,removal:end);

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
    disp('tapering segmentation processed')
    
    %% identify start and terminal node
    s_all_term_idx = find([snode.ep]);
    sK = ceil([snode(s_all_term_idx).comz]); % vox indicies of terminal points
    %[sJ,sI,sK] = ind2sub(size(s_skel),sterminal_vox);

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
            [~, max_idx] = max(raw_pks(min_idx-search_area:min_idx+search_area));
            idx = min_idx - (search_area+1 - max_idx);
        catch
             idx = min_idx;
        end
        raw_pks(idx) = 0; % prevents same peak being chosen twice.
        new_bifurcation_idx(i) = raw_locs(idx);
        % consider something here incase finding peak fails
    end
    disp('Bifurcation points identified')

    %% Remove Bifurcation Points
    % most computationally intensive.
    % add carina to be excluded.
    new_bifurcation_idx = [1 new_bifurcation_idx];

    %% save result
    disp('saving results')
    N1_cleaned_data(j).bi_loc = new_bifurcation_idx;
    catch
    end
end
