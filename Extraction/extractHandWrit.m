function [HandWrit,Skew,Frame] = extractHandWrit(RGB,SqrSz)
%EXTRACT return RGB image of handwritten
% Input:
%       RGB: RGB image
%       SqrSz: size of square indicator
% Output:
%       HandWrit: RGB image of handwritten
%       Skew: skew of image
%       Frame: frame to cut handwritten

% convert to black and white temporary
Temp = im2bw(RGB);
Temp = ~Temp;% black is true

SqrCC = findSqr(Temp,SqrSz);
% remove squares 
for i=1:SqrCC.NumObjects
    R = RGB(:,:,1);
    G = RGB(:,:,2);
    B = RGB(:,:,3);
    R(SqrCC.PixelIdxList{i}) = intmax(class(RGB));
    G(SqrCC.PixelIdxList{i}) = intmax(class(RGB));
    B(SqrCC.PixelIdxList{i}) = intmax(class(RGB));
    RGB = cat(3,R,G,B);
end

Mask = ROIPolygon(SqrCC);

Ind = find(Mask);
% remove connected component on border of Mask
HandWrit = ind2img(Ind,RGB);
Skew = SkewAngle(SqrCC);
Frame = bwperim(Mask);
FrameIndex = find(Frame);
Frame = ind2img(FrameIndex,Frame);

end

