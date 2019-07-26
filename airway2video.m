% Script to visualise taper results.

%% load
%load('P_N1_raw_image_CT_nifti_A_22195509.mat')

%% extract and reassign variables
cross_images = [sturct_to_save.tapering_raw_image.resample_image];
info = cell2mat(sturct_to_save.tapering_raw_image.area_results.elliptical_info);
wall_info = cell2mat(sturct_to_save.tapering_raw_image.area_results.elliptical_info_wall);
lumen_x = {info(:).x_stop}';
lumen_y = {info(:).y_stop}';
wall_x = {wall_info(:).x_stop}';
wall_y = {wall_info(:).y_stop}';
lumen_area = [sturct_to_save.tapering_raw_image.area_results.phyiscal_area];
wall_area = [sturct_to_save.tapering_raw_image.area_results.phyiscal_area_wall];
arc_length = [sturct_to_save.tapering_raw_image.arclegth];
lumen_ellipses = {info(:).elipllical_info}';
wall_ellipses = {wall_info(:).elipllical_info}';

%% set up video object
video = VideoSetup('airwayN1_22195509_peak_anom_med');
video.open();

%% create video frames
h = figure;
for i = 1:size(cross_images,3)
    % display image
    imagesc(cross_images(:,:,i))
    hold on
    colormap gray
    
    % plot ray cast results
    plot(lumen_x{i, 1},lumen_y{i, 1},'r.')
    plot(wall_x{i, 1},wall_y{i, 1},'c.')
    
    % plot ellipse fitting
    ellipse(lumen_ellipses{i,1}(3),lumen_ellipses{i,1}(4),...
        lumen_ellipses{i,1}(5),lumen_ellipses{i,1}(1),...
        lumen_ellipses{i,1}(2),'m')
    
    ellipse(wall_ellipses{i,1}(3),wall_ellipses{i,1}(4),...
        wall_ellipses{i,1}(5),wall_ellipses{i,1}(1),...
        wall_ellipses{i,1}(2),'b')

    
    % display area measurements
    dim = [.15 .85 .24 .05];
    %a = annotation('textbox',dim,'String',str,'FitBoxToText','on','BackgroundColor','y');
    a = rectangle('Position',[0,0,133,10],'FaceColor','y','LineWidth',2);
    ax = gca;
    txt = text(ax, 1,5,sprintf('Arc Length = %4.1f mm; Lumen area = %4.2f mm^2; Wall area = %4.2f mm^2; %3.0i of %3.0i', ...
        arc_length(i), lumen_area(i),wall_area(i), i, size(cross_images,3)));
    
    
    drawnow
    
    VideoAddFrame(video, h);

    delete(a)
    delete(txt)
    
end

%% save
video.close();

