% By Ashkan Pakzad on 30/07/2019 

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
