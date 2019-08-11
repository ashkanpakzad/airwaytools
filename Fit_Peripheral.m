% Fitting peripheral taper rate.
% by Ashkan Pakzad on 10th Aug 2019

prefix = 'P_S1_raw_image_CT_nifti_A_';
% create list of airway ids
airway_number = [S1_cleaned_data(:).airway_number];

%% run loop to determine bifurcation data points

for j = 1:length(airway_number)
    disp(['airway ', num2str(j) ' of ' , num2str(length(airway_number))])
    try % incase airway missing
        %% load
        filename = [prefix,num2str(airway_number(j)),'.mat'];

        % extract info from airway struct.
        arc_length = S1_cleaned_data(j).arc_length;
        lumen_area = S1_cleaned_data(j).lumen_area;
        wall_area = S1_cleaned_data(j).wall_area;
        bifurcation_idx = S1_cleaned_data(j).bi_loc;
        logic_include = S1_cleaned_data(j).nonbi_include;

        %% peripheral taper fitting
        % % fitting from the xth bifurcation (not including carina).
        bifurcation_fit = 5;

        [p_lumen, p_wall, fit_arc_length] = PeripheralTaper(bifurcation_fit,bifurcation_idx, logic_include, arc_length, lumen_area,wall_area);
        assert((p_lumen(1) ~= 0) && (p_wall(1) ~= 0))
        line_lumen = p_lumen(1)*fit_arc_length+p_lumen(2);
        line_wall = p_wall(1)*fit_arc_length+p_wall(2);

        
        %% fit taper rate
        [ log_taper_rate_lumen, log_taper_rate_wall ] =  ...
        Log_taper_rate_of_ordered_array(arc_length(logic_include),lumen_area(logic_include),wall_area(logic_include));

        %% save
        S1_cleaned_data(j).lumen_t5 = p_lumen;
        S1_cleaned_data(j).wall_t5 = p_wall;
        S1_cleaned_data(j).fit_arc_length_t5 = fit_arc_length;
        
        S1_cleaned_data(j).lumen_log_taper_rate = log_taper_rate_lumen;
        S1_cleaned_data(j).wall_log_taper_rate = log_taper_rate_wall;    


    catch
         S1_cleaned_data(j).lumen_t5 = [];
         S1_cleaned_data(j).wall_t5 = [];
         S1_cleaned_data(j).fit_arc_length_t5 = [];
     end
    if S1_cleaned_data(j).lumen_log_taper_rate == 0 | S1_cleaned_data(j).wall_log_taper_rate == 0
        S1_cleaned_data(j).lumen_log_taper_rate = [];
        S1_cleaned_data(j).wall_log_taper_rate = [];
    end
end