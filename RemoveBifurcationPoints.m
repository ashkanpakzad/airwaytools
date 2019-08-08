function [non_bi_points, new_bifurcation_idx, false_counter] = RemoveBifurcationPoints(removal, search_area, s_raw, lumen_area)
%% INPUT
%    removal = number of initial values to ignore in lumen area in order to
%    remove carina from skeletisation of the tapering segmentation.
%    recommend = 50.
%    search_area = the half-width of peaks to consider when considering
%    local area peaks as bifurcations. recommend = 4.
%    s_raw = resampled segmentation from airway tapertool.
%    lumen_area = array of lumen area from airway tapertool.

[sG, snode, slink] = TaperingSkeleton(s_raw,removal);
[bifurcation_idx,false_counter] = FindBifurcationsSkel(sG, snode, slink);
[new_bifurcation_idx] = FindBifurcationsPeak(lumen_area, bifurcation_idx, search_area);
non_bi_points = RemoveBifurcations(new_bifurcation_idx);

%%
function [sG, snode, slink] = TaperingSkeleton(s_raw,removal)
% Preprocessing of indiv. airway segmentation
s_raw = s_raw(:,:,removal:end);

% force binary image
sraw_binary = s_raw > 0.7;

% smoothen binary image
N = 3;
kernel = ones(N, N, N) / N^3;
blurryImage = convn(double(sraw_binary), kernel, 'same');
s_binary = blurryImage > 0.5;

% fill holes
s_fill = imfill(s_binary, 'holes');

% remove smaller structures.
s = zeros(size(s_fill));
CC = bwconncomp(s_fill);
numPixels = cellfun(@numel,CC.PixelIdxList); % apply func to each cell
[~,idx] = max(numPixels);
s(CC.PixelIdxList{idx}) = 1;
s = logical(s);

% skeletonise and convert single airway to graph
s_skel = bwskel(s);

[sA,snode,slink] = Skel2Graph3D(s_skel,0); % dependent function
sG = graph(sA);
end

%%

function [bifurcation_idx,false_counter] = FindBifurcationsSkel(sG, snode, slink)
% identify start and terminal node
s_all_term_idx = find([snode.ep]);
sK = ceil([snode(s_all_term_idx).comz]); % vox indicies of terminal points
%[sJ,sI,sK] = ind2sub(size(s_skel),sterminal_vox);

% first z node = start
[~,start_zidx] = min(sK');
start_idx = s_all_term_idx(start_zidx);

% last z node = airway end
[~,stop_zidx] = max(sK');
stop_idx = s_all_term_idx(stop_zidx);

sP = shortestpath(sG,start_idx,stop_idx); % path for airway

% remove false bifurcations by analysing skeleton
% Func. to find plausible bifurcations
[true_bifurcations, ~] = Bifurcation_Plausibility(sP, snode, slink);
% extract plausible bifurcations
sP_bi = sP(2:end-1); % removes start and end nodes
sP_plaus = sP_bi(true_bifurcations);
% re-add start and end nodes
sP_plaus = [sP(1), sP_plaus, sP(end)];

% counter for false bifurcations
false_counter = sum(~true_bifurcations);

% identify arc_length distances for path bifurcation nodes
path_zind = round([snode(sP_plaus).comz]);
%[pathJ,pathI,pathK] = ind2sub(size(s_skel),path_ind);
path_zind = path_zind(2:end-1); % remove start and end
bifurcation_idx = path_zind+removal;

    function [true_bifurcations, link_len] = Bifurcation_Plausibility(Path,node, link)
    % NB: must contain first and last points (which are terminals)

    % Path nodes from carina to terminal point
    % node and link from skel2graph

    % output = vector of size path-2. lengths of links to alternate path.
    % long length = true bifurcation (generally).
    % By Ashkan Pakzad 27th July 2019

    link_len = zeros(length(Path)-2,1);
    for i = 2:length(Path)-1
        connected_nodes = node(Path(i)).conn;
        connected_links = node(Path(i)).links;
        off_branch_idx = find(connected_nodes ~= Path(i-1) & ...
            connected_nodes ~= Path(i+1));
        off_branch_link = connected_links(off_branch_idx);
        link_len(i-1) = numel([link(off_branch_link).point]);
    end

    %[~,L,~,~] = isoutlier(link_len,'mean');
    true_bifurcations = link_len > 3;
    end
end
%%

function [new_bifurcation_idx] = FindBifurcationsPeak(lumen_area, bifurcation_idx, search_area)
% identify bifurcations with local maxima peaks raw
[raw_pks,raw_locs] = findpeaks(lumen_area);
new_bifurcation_idx = zeros(size(bifurcation_idx));
for i = 1:length(bifurcation_idx)
    % identify nearest peak
    [~, min_idx] = min(abs(raw_locs - bifurcation_idx(i)));
    try
        % consider if nearer peaks are more likely to rep. bifurcation
        [~, max_idx] = max(raw_pks(min_idx-search_area:min_idx+search_area));
        idx = min_idx - (search_area+1 - max_idx);
    catch
         idx = min_idx;
    end
    raw_pks(idx) = 0; % prevents same peak being chosen twice.
    new_bifurcation_idx(i) = raw_locs(idx);
    % consider something here incase finding peak fails
end
end

function non_bi_points = RemoveBifurcations(new_bifurcation_idx)
%Remove Bifurcation Points
% most computationally intensive.
% add carina to be excluded.
new_bifurcation_idx = [1 new_bifurcation_idx];
%indices for exclusion between lower and upper point at each bifurcation.
lower_exclude = zeros(size(new_bifurcation_idx));
upper_exclude = zeros(size(new_bifurcation_idx));

for i = 1:length(new_bifurcation_idx)

    pck_idx = new_bifurcation_idx(i);
    try
        if i == 1
            data_hwidth = 100;% initial peak width for carina
        else
            data_hwidth = 50;% initial peak width
        end
        % initiate loop variables
        gof_n = 0;
        gof_np1 = 0.0001;
        int = 1; % inteval width reduction
        lower_width = pck_idx-data_hwidth;
        upper_width = pck_idx+data_hwidth;
        % set precautions incase 1 < data width > last idx.
        if lower_width < 1
            lower_width = 1;
        else
        end
        if upper_width > length(lumen_area)
            upper_width = length(lumen_area);
        else
        end
        x_data = mat2gray(lumen_area(lower_width:upper_width));
        [~, x_data_pck_idx] = max(x_data);

         % reduce data window if bifurcation is not the peak.
        if i ~= 1
            while x_data_pck_idx ~= data_hwidth+1
                [~, x_data_pck_idx] = max(x_data);
                data_hwidth = data_hwidth - int;
                if data_hwidth == 3
                    break
                end
            end
        else
        end

        while gof_n < gof_np1 || gof_n < 0.2
            try % save prev. result for display.
                %n_fitresult = fitresult; 
                n_data_hwidth = data_hwidth;
                %n_x_data = x_data;
                %n_vals = vals;
                catch
            end

            % normalise data window
            if lower_width < 1
                lower_width = 1;
            else
            end
            if upper_width > length(lumen_area)
                upper_width = length(lumen_area);
            else
            end
            x_data = mat2gray(lumen_area(lower_width:upper_width));

            % equally spaced values
            vals = [1:length(x_data)]';

            % fit gaussian with fixed mean at centre
            options = fitoptions('gauss1');
            if i == 1 % mean is 1 for carina
                options.Lower = [-inf 1 0];
                options.Upper = [inf 1 inf];
            else
                options.Lower = [-inf data_hwidth+1 0];
                options.Upper = [inf data_hwidth+1 inf];
            end
            [fitresult, gof] = fit(vals, x_data, 'gauss1', options);
            %[~, gof] = fit(vals, x_data, 'gauss1', options);


            % assign values for next loop
            gof_n = gof_np1; 
            gof_np1 = gof.rsquare;
            % reduce peak width 
            data_hwidth = data_hwidth - int;
            lower_width = pck_idx-data_hwidth;
            upper_width = pck_idx+data_hwidth;
        end
    catch
        n_data_hwidth = 3; % if fail to fit, remove 7 points around bifurcation
    end
    if i == 1
        lower_exclude(i) = 1;
    else
        lower_exclude(i) = pck_idx-n_data_hwidth;
    end
    upper_exclude(i) = pck_idx+n_data_hwidth;
end
% extract all points not excluded
non_bi_points = true(size(lumen_area));
for i = 1:length(lower_exclude)
    non_bi_points(lower_exclude(i):upper_exclude(i)) = 0;
end
end
end