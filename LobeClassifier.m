% By Ashkan Pakzad 30/07/19: lobe classifier of airway terminal points.

function [airway_number, classification] = LobeClassifier(airway_seg_name, lobe_seg_name)
% classification labels
% Right Superior    - 1
% Right Middle      - 2
% Right Inferior    - 3
% Left Superior     - 4
% Left Inferior     - 5

L = niftiread(lobe_seg_name);

% fill lobe holes
F = zeros(size(L));
for i = 1:size(L,3)
    F(:,:,i) = imfill(L(:,:,i));
end

DistalPoints(airway_seg_name,'distal',0)
D = niftiread('distal');
airway_number = find(D);

classification = F(airway_number);
end