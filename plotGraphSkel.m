function plotGraphSkel(skel, node, type)

% Function by Ashkan Pakzad on 27th July 2019 for use after Skel2Graph.

if nargin < 3
  inputtype = 'full';
else
  inputtype = type;
end

if strcmp(inputtype, 'full') == 1
    X = [node.comy];
    Y = [node.comx];
    Z = [node.comz];
    nums = string(1:length(X));
    
elseif strcmp(inputtype, 'terminals') == 1
    term_idx = find([node.ep]);
    X = [node(term_idx).comy];
    Y = [node(term_idx).comx];
    Z = [node(term_idx).comz];
    nums = string(node(term_idx));
    
else
    error('Type must be ''full'' or ''terminals''')
end

isosurface(skel);
hold on
plot3(X,Y,Z, 'r.', 'MarkerSize', 15);
text(X+1,Y+1,Z+1, nums)