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

%% display box plot per lobe N1 only
Lobe_idxN1= cell(1,5);
Lobe_taperN1 = cell(1,5);
box_dataN1 = [];
grpN1 = [];
label = ["RS","RM","RI","LS","LI"];

for i = 1:5 % for all lobes
    Lobe_idxN1{1, i} = find([N1_cleaned_data(:).class] == i);
    Lobe_taperN1{1, i} = [N1_cleaned_data(Lobe_idxN1{1, i}).lumen_log_taper_rate];
    box_dataN1 = [box_dataN1; Lobe_taperN1{1, i}'];
    grpN1 = [grpN1; repmat(label(i),length(Lobe_taperN1{1, i}),1)];
end

% ["RS","RM","RI","LS","LI"]

figure
%boxplot(box_data,grp)
boxplot2(box_data,grp)
title('N1 Lumen Log Taper rate per lobe')
xlabel('Lobe')
ylabel('Log Taper rate (mm^-1)') % taking log makes it dimensionless.

%% display box plot per lobe M1

Lobe_idxM1 = cell(1,5);
Lobe_taperM1 = cell(1,5);
box_dataM1 = [];
grpM1 = [];
label = ["RS","RM","RI","LS","LI"];
for i = 1:5 % for all lobes
    Lobe_idxM1{1, i} = find([M1_cleaned_data(:).class] == i);
    Lobe_taperM1{1, i} = [M1_cleaned_data(Lobe_idxM1{1, i}).lumen_log_taper_rate];
    box_dataM1 = [box_dataM1; Lobe_taperM1{1, i}'];
    grpM1 = [grpM1; repmat(label(i),length(Lobe_taperM1{1, i}),1)];
end

% ["RS","RM","RI","LS","LI"]

figure
boxplot(box_dataM1,grpM1)
title('M1 Lumen Log Taper rate per lobe')
xlabel('Lobe')
ylabel('Log Taper rate (mm^-1)') % taking log makes it dimensionless.

%% cluster boxplot N1 M1
x = 1:5; % 5 clusters
maxval = 84; % max number of values in one boxplot
y = nan(5, 2, maxval); % numberofclusters x plot in each cluster x data

for i = 1:5
    y(i, 1, 1:length(Lobe_taperN1{1,i})) = Lobe_taperN1{1,i};
    y(i, 2, 1:length(Lobe_taperM1{1,i})) = Lobe_taperM1{1,i};
end

% Plot boxplots
figure
h = boxplot2(y,x);
ax = gca;
ax.XTick = x;
ax.XTickLabel = ["RS","RM","RI","LS","LI"];
title('Lumen Log Taper rate per lobe')
xlabel('Lobe')
ylabel('Log Taper rate (mm^{-1})') % taking log makes it dimensionless.

cmap = get(0, 'defaultaxescolororder');
for ii = 1:2
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end

%legend(findall(ax,'Tag','Box'),{'M1', 'N1'});
legend(ax,{'N1', 'M1'});


%% stat tests
p = zeros(5,1);
for i = 1:5
    p(i) = ranksum(Lobe_taperN1{1,i},Lobe_taperM1{1,i});
end
disp(label)
disp(p')