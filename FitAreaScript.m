% By Ashkan Pakzad on 10th August 2019
% Find lumen edge and wall attenuation peak, fits ellipses to both and
% computes area.

% intended for use with S1, block 2 airways only.
%% Load CT image and extract parameters
[raw_header,~,~] = Complete_image_load_nii('S1.nii');
plane_input_sturct = Return_airway_interpol_para(raw_header);

%% Load all airway names
% Specify the folder where the files live.
myFolder = '../blocks2';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

%% Process airways one by one
for i = 1 : length(theFiles)
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    disp(['Airway ', num2str(i),' of ',num2str(length(theFiles))])
    load(fullFileName);

    area_results = Computing_lumen_FWHM_SL_with_safety_catch(...
            sturct_to_save.tapering_raw_image.resample_image,...
            plane_input_sturct,...
            sturct_to_save.tapering_seg_image.resample_image);

    %compute the tapering measurment
    [area_results.tapering, area_results.tapering_wall] = ...
            Log_taper_rate_of_ordered_array(...
            sturct_to_save.tapering_raw_image.arclegth,...
            area_results.phyiscal_area, area_results.phyiscal_area_wall);

    %% save to same convention
    sturct_to_save.tapering_raw_image.area_results = area_results;
    savename = baseFileName(1:end-4);
    save(savename,'sturct_to_save') % save without .mat ext.
    clearvars sturct_to_save % clear to prep for next airway
end