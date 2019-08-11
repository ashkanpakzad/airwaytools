function [box_data, grp, Lobe_taper] = ResultsFigs(cleaned_data, measuretype)
% By Ashkan Pakzad on  11th August 2019
% Groups individual subject data to produce box plot.

Lobe_idx= cell(1,5);
Lobe_taper = cell(1,5);
box_data = [];
grp = [];
label = ["RU","RM","RL","LU","LL"];

if strcmp(measuretype, 'peripheral') == 1
    for i = 1:5
        Lobe_idx{1, i} = find([cleaned_data(:).class] == i);
        p_store = [cleaned_data(Lobe_idx{1, i}).lumen_t5];
        Lobe_taper{1, i} = p_store(1:2:end);
        box_data = [box_data; Lobe_taper{1, i}'];
        grp = [grp; repmat(label(i),length(Lobe_taper{1, i}),1)];
    end

elseif strcmp(measuretype, 'log') == 1
    for i = 1:5
        Lobe_idx{1, i} = find([cleaned_data(:).class] == i);
        Lobe_taper{1, i} = [cleaned_data(Lobe_idx{1, i}).lumen_log_taper_rate];
        box_data = [box_data; Lobe_taper{1, i}'];
        grp = [grp; repmat(label(i),length(Lobe_taper{1, i}),1)];
    end
else
    error('Type must be ''peripheral'' or ''log'' of taper.')
end
end
