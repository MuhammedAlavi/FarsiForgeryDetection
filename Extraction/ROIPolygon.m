function Mask = ROIPolygon(SqrCC)
% ROIPOLYGON Region Of Interest according to 4 square indicator
% Input:
%       SqrCC: connected component structure of 4 square indicator
% Output:
%       Mask: binary mask of region of interest
s = cellfun(@sum,SqrCC.PixelIdxList,'UniformOutput',false);
s = cell2mat(s);
[~,I] = sort(s,'ascend');

% find center of each indicator square
center = [0,0];
Temp = zeros(SqrCC.ImageSize);
for i=1:SqrCC.NumObjects
    Temp(SqrCC.PixelIdxList{i}) = 1;
    c = regionprops(Temp,'centroid');
    center(i,:) = c.Centroid(:);
    Temp = zeros(size(Temp));
end

Temp = zeros(size(Temp));
Polygon = [center(I(1),:);center(I(3),:);center(I(4),:);...
           center(I(2),:);center(I(1),:)];
Mask = roipoly(Temp,Polygon(:,1),Polygon(:,2));
end