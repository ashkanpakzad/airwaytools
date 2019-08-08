% Script to process S1 airway tapering.
% By Ashkan Pakzad on 8th July 2019

%% input var
prefix = 'P_S1_raw_image_CT_nifti_A_';
% remove bronchi at carina
removal = 50;
search_area = 4; % for local maxima

%% get airway numbers, classification & initiate results
[airway_number, classification] = LobeClassifier('S1_raw_image_lumen_seg.nii', 'S1_lobe_seg.nii.gz');

S1_cleaned_data = struct('airway_number', [],...
    'class',[],'arc_length',{}, 'lumen_area',{}, 'wall_area',{}, ...
    'nonbi_include',{},'lumen_log_taper_rate',[], 'wall_log_taper_rate', []);

%% run loop to determine bifurcation data points
for j = 1:length(airway_number)
    disp(['airway ', num2str(j) ' of ' , num2str(length(airway_number))])
    
    %% load
    filename = [prefix,num2str(airway_number(j)),'.mat'];
    load(filename);
    s_raw = sturct_to_save.tapering_seg_image.resample_image;
    
    % extract info from airway struct.
    arc_length = sturct_to_save.tapering_raw_image.arclegth;
    area_info = sturct_to_save.tapering_raw_image.area_results;
    lumen_area = [area_info.phyiscal_area];
    wall_area = [area_info.phyiscal_area_wall];
    
    %% remove bifurcations
    [logic_include, false_counter] = RemoveBifurcationPoints(removal, search_area, s_raw, lumen_area);
    
%     bi_removed_arcl= arc_length(logic_include);
%     bi_removed_lumen_area = lumen_area(logic_include);
%     bi_removed_wall_area = wall_area(logic_include);

    %% fit petaperline
%     [ log_taper_rate_lumen, log_taper_rate_wall ] =  ...
%     Log_taper_rate_of_ordered_array(bi_removed_arcl,bi_removed_lumen_area,bi_removed_wall_area);

    %% save result
    disp('saving results')
    S1_cleaned_data(j).airway_number = airway_number(j);
    S1_cleaned_data(j).class = classification(j);
    S1_cleaned_data(j).arc_length = arc_length;
    S1_cleaned_data(j).lumen_area = lumen_area;
    S1_cleaned_data(j).wall_area = wall_area;
    S1_cleaned_data(j).nonbi_include = logic_include;
    S1_cleaned_data(j).lumen_log_taper_rate = log_taper_rate_lumen;
    S1_cleaned_data(j).wall_log_taper_rate = log_taper_rate_wall;
    
    %% save result
    save('S1_cleaned_data', 'S1_cleaned_data')
end

