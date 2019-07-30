function plotGraphSkel(skel, node)

% Function by Ashkan Pakzad on 27th July 2019 for use after Skel2Graph.

isosurface(skel);
hold on
plot3([node.comy],[node.comx], [node.comz], 'r.', 'MarkerSize', 15);
text([node.comy]+1,[node.comx]+1, [node.comz]+1, string(1:length(node)))