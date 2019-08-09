% testing ground for removing bifurcations and plotting results.
% by Ashkan Pakzad on 8th Aug 2019

%% load
filename = 'P_S1_raw_image_CT_nifti_A_52042583.mat';
load(filename);
s_raw = sturct_to_save.tapering_seg_image.resample_image;

% extract info from airway struct.
arc_length = sturct_to_save.tapering_raw_image.arclegth;
area_info = sturct_to_save.tapering_raw_image.area_results;
lumen_area = [area_info.phyiscal_area];
wall_area = [area_info.phyiscal_area_wall];

%% remove bifurcations
removal = 50;
search_area = 7;
[logic_include, bifurcation_idx, false_counter] = RemoveBifurcationPoints(removal, search_area, s_raw, lumen_area);

bi_removed_arcl = arc_length(logic_include);
bi_removed_lumen_area = lumen_area(logic_include);
bi_removed_wall_area = wall_area(logic_include);

%% peripheral taper fitting
% % fitting from the xth bifurcation (not including carina).
bifurcation_fit = 5;

[p_lumen, p_wall, fit_arc_length] = PeripheralTaper(bifurcation_fit,bifurcation_idx, logic_include, arc_length, lumen_area,wall_area);
line_lumen = p_lumen(1)*fit_arc_length+p_lumen(2);
line_wall = p_wall(1)*fit_arc_length+p_wall(2);

%% visualise
C = linspecer(5,'sequential'); % distinct colour creator by Jonathan Lansey

figure
plot(arc_length, lumen_area,'Color',C(1,:),'Marker','.','LineStyle','none');
hold on
plot(bi_removed_arcl, bi_removed_lumen_area, 'Marker','.','Color',C(2,:),'LineStyle','none');
plot(arc_length(bifurcation_idx), lumen_area(bifurcation_idx), 'Color',C(3,:),'Marker','^','LineStyle','none','MarkerSize',10,'MarkerFaceColor',C(4,:))
plot(fit_arc_length, line_lumen,'Color',C(5,:),'LineWidth',2)
title(['Lumen area - Arc Length \tau_{', num2str(bifurcation_fit),'}'])
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Measurements','Non-Bifurcation Measurements', 'Bifurcation Peak', 'Peripheral tapering fit')
text(fit_arc_length(1)+20, line_lumen(1)+50, ['\tau_{', num2str(bifurcation_fit),'} = ',num2str(p_lumen(1))],'Color',C(5,:),'FontSize',14)
axis tight

%% visualise wall

figure
plot(arc_length, wall_area,'Color',C(1,:),'Marker','.','LineStyle','none');
hold on
plot(bi_removed_arcl, bi_removed_wall_area, 'Marker','.','Color',C(2,:),'LineStyle','none');
plot(arc_length(bifurcation_idx), wall_area(bifurcation_idx), 'Color',C(3,:),'Marker','^','LineStyle','none','MarkerSize',10,'MarkerFaceColor',C(4,:))
plot(fit_arc_length, line_wall,'Color',C(5,:),'LineWidth',2)
title(['Wall area - Arc Length \tau_{', num2str(bifurcation_fit),'}'])
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Measurements','Non-Bifurcation Measurements', 'Bifurcation Peak', 'Peripheral tapering fit')
text(fit_arc_length(1)+20, line_wall(1)+50, ['\tau_{', num2str(bifurcation_fit),'} = ',num2str(p_wall(1))],'Color',C(5,:),'FontSize',14)
axis tight
