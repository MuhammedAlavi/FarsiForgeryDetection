function d = SkewAngle(cc)
%SKEWANGLE determine skew angle of image
% Input:
%       cc: connected component of 4 square image
% output: 
%       d: degree of skew

% add a new field to cc
cc.SumOfRows = zeros(1,cc.NumObjects);

% convert to subsctip
for i=1:cc.NumObjects
    [Row,Col] = ind2sub(cc.ImageSize,cc.PixelIdxList{i});
    cc.PixelIdxList{i} = [Row,Col];
    cc.SumOfRows(i) = sum(Row);
end

% find two upper square(two square with minimum sum of rows)
[~,I] = sort(cc.SumOfRows,'ascend');
LU = cc.PixelIdxList{I(1)};% left upper
RU = cc.PixelIdxList{I(2)};% right upper

% difference of rows and columns for first point
Diff = LU(1,:) - RU(1,:);

% degree of skew
d = atan2d(Diff(1),Diff(2));

if d < -90
    d = 180 + d;
end
end

