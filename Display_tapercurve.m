% By Ashkan Pakzad on 30/07/2019 

%% display lumen taper curves
j = 1;
arc = N1_cleaned_data(j).arc_length;
area = N1_cleaned_data(j).lumen_area;
include = N1_cleaned_data(j).nonbi_include;

newarc = arc(include);
newarea = area(include);

area_log = log(newarea);

% fit w/o bifurcation points
p_coeff = polyfit(newarc,area_log,1);
plot_data = exp(p_coeff(1)*arc + p_coeff(2));

% fit w/ bifurcation points
p_coeff_all = polyfit(arc,log(area),1);
plot_data_all = exp(p_coeff_all(1)*arc + p_coeff_all(2));

figure
plot(arc, area, '.');
hold on
plot(newarc, newarea, '.');
plot(arc, plot_data_all, '-');
plot(arc, plot_data, '-');
title('Lumen area - Arc Length')
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Points','Non-Bifurcation points','Taper Curve w/ BP', 'Taper Curve w/o BP')

%% display wall taper curves
j = 1;
arc = N1_cleaned_data(j).arc_length;
area = N1_cleaned_data(j).wall_area;
include = N1_cleaned_data(j).nonbi_include;

newarc = arc(include);
newarea = area(include);

area_log = log(newarea);

% fit w/o bifurcation points
p_coeff = polyfit(newarc,area_log,1);
plot_data = exp(p_coeff(1)*arc + p_coeff(2));

% fit w/ bifurcation points
p_coeff_all = polyfit(arc,log(area),1);
plot_data_all = exp(p_coeff_all(1)*arc + p_coeff_all(2));

figure
plot(arc, area, '.');
hold on
plot(newarc, newarea, '.');
plot(arc, plot_data_all, '-');
plot(arc, plot_data, '-');
title('Wall area - Arc Length')
xlabel('Arc Length (mm)')
ylabel('Area (mm^2)')
legend('Bifurcation Points','Non-Bifurcation points','Taper Curve w/ BP', 'Taper Curve w/o BP')

%% display box plot per lobe
Lobe_idx = cell(1,5);
Lobe_taper = cell(1,5);
box_data = [];
grp = [];
label = ["RS","RM","RI","LS","LI"];

for i = 1:5 % for all lobes
    Lobe_idx{1, i} = find([N1_cleaned_data(:).class] == i);
    Lobe_taper{1, i} = [N1_cleaned_data(Lobe_idx{1, i}).lumen_log_taper_rate];
    box_data = [box_data; Lobe_taper{1, i}'];
    grp = [grp; repmat(label(i),length(Lobe_taper{1, i}),1)];
end

% ["RS","RM","RI","LS","LI"]

figure
boxplot(box_data,grp)
title('Lumen Log Taper rate per lobe')
xlabel('Lobe')
ylabel('Log Taper rate (mm^-1)') % taking log makes it dimensionless.