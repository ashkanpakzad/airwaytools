% By Ashkan Pakzad on 30th July 2019
%% display box plot per lobe N1
[box_dataN1, grpN1, ~] = ResultsFigs(N1_cleaned_data, 'peripheral');

figure
boxplot(box_dataN1,grpN1)
title('N1 Lumen Peripheral Taper rate per lobe')
xlabel('Lobe')
ylabel('Taper rate (mm)') % taking log makes it dimensionless.

%% display box plot per lobe M1 Peripheral taper
[box_dataM1, grpM1, ~] = ResultsFigs(M1_cleaned_data, 'peripheral');

figure
boxplot(box_dataM1,grpM1)
title('M1 Lumen Peripheral Taper rate per lobe')
xlabel('Lobe')
ylabel('Taper rate (mm)') % taking log makes it dimensionless.

%% display box plot per lobe S1 Peripheral taper
[box_dataS1, grpS1, ~] = ResultsFigs(M1_cleaned_data, 'peripheral');

figure
boxplot(box_dataS1,grpS1)
title('S1 Lumen Peripheral Taper rate per lobe')
xlabel('Lobe')
ylabel('Taper rate (mm)') % taking log makes it dimensionless.

%% Load all vars.
load('N1_cleaned_data_v4.mat')
load('M1_cleaned_data_v4.mat')
load('S1_cleaned_data_v4.mat')

%% cluster boxplot N1 M1 S1
x = 1:5; % 5 clusters
maxval = 135; % max number of values in one boxplot
y = nan(5, 3, maxval); % numberofclusters x plot in each cluster x data

type = 'log';
[~, ~, Lobe_taperN1] = ResultsFigs(N1_cleaned_data, type);
[~, ~, Lobe_taperM1] = ResultsFigs(M1_cleaned_data, type);
[~, ~, Lobe_taperS1] = ResultsFigs(S1_cleaned_data, type);

for i = 1:5
    y(i, 1, 1:length(Lobe_taperN1{1,i})) = Lobe_taperN1{1,i};
    y(i, 2, 1:length(Lobe_taperM1{1,i})) = Lobe_taperM1{1,i};
    y(i, 3, 1:length(Lobe_taperS1{1,i})) = Lobe_taperS1{1,i};
end

% Plot boxplots
figure
h = boxplot2(y,x);
ax = gca;
ax.XTick = x;
ax.XTickLabel = ["RU","RM","RL","LU","LL"];
title('Lumen Log Taper rate per lobe')
xlabel('Lobe')
ylabel('Log Taper rate (mm^{-1})') % taking log makes it dimensionless.

cmap = get(0, 'defaultaxescolororder');
for ii = 1:3
    structfun(@(x) set(x(ii,:), 'color', cmap(ii,:), ...
        'markeredgecolor', cmap(ii,:)), h);
end

%legend(findall(ax,'Tag','Box'),{'M1', 'N1'});
legend(ax,{'N1', 'M1', 'S1'});


%% stat tests
p = zeros(5,1);
for i = 1:5
    p(i) = ranksum(Lobe_taperN1{1,i},Lobe_taperM1{1,i});
end
disp(label)
disp(p')