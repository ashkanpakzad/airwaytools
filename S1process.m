% Script to save N1 bi points.
% By Ashkan Pakzad on 10th August 2019

%% input var
data_folder = fullfile('..','..','results','S1');
prefix = 'P_S1_raw_image_CT_nifti_A_';
% remove bronchi at carina
removal = 50;
search_area = 7; % for local maxima

%% get airway numbers, classification & initiate results
[airway_number, classification] = LobeClassifier('S1_clumen.nii', 'S1_lobe_seg_need_check.nii');

S1_cleaned_data = struct('airway_number', [],...
    'class',[],'arc_length',{}, 'lumen_area',{}, 'wall_area',{}, ...
    'nonbi_include',{},'bi_loc',{}, 'no_false_bifurcations', []);

%% run loop to determine bifurcation data points
for j = 1:length(airway_number)
    disp(['airway ', num2str(j) ' of ' , num2str(length(airway_number))])
    try % incase airway missing
    %% load
    filename = [prefix,num2str(airway_number(j)),'.mat'];
    load(fullfile(data_folder,filename));
    s_raw = sturct_to_save.tapering_seg_image.resample_image;
    
    % extract info from airway struct.
    arc_length = sturct_to_save.tapering_raw_image.arclegth;
    area_info = sturct_to_save.tapering_raw_image.area_results;
    lumen_area = [area_info.phyiscal_area];
    wall_area = [area_info.phyiscal_area_wall];

    % remove bifurcations
    [logic_include, bifurcation_idx, false_counter] = RemoveBifurcationPoints(removal, search_area, s_raw, lumen_area);

    bi_removed_arcl = arc_length(logic_include);
    bi_removed_lumen_area = lumen_area(logic_include);
    bi_removed_wall_area = wall_area(logic_include);

    % save result
    S1_cleaned_data(j).airway_number = airway_number(j);
    S1_cleaned_data(j).class = classification(j);
    S1_cleaned_data(j).arc_length = arc_length;
    S1_cleaned_data(j).lumen_area = lumen_area;
    S1_cleaned_data(j).wall_area = wall_area;
    S1_cleaned_data(j).nonbi_include = logic_include;
    S1_cleaned_data(j).bi_loc = bifurcation_idx;
    S1_cleaned_data(j).no_false_bifurcations = false_counter;

    catch
        S1_cleaned_data(j).airway_number = airway_number(j);
        warning(['Airway ', num2str(airway_number(j)), ' is missing.'])
    end
    %% save result
    save('S1_cleaned_data_v2', 'S1_cleaned_data')
end
