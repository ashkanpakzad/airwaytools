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

[~,L,~,~] = isoutlier(link_len);
true_bifurcations = link_len > L;
end
    