function ClearCC = removeOnBorder(BW,CC,Frame)
%REMOVEONBORDER remove connected components on border 
% Input: 
%       BW: black and white image
%       CC: connected component structure
%       Frame: frame when crop handwritten
% Output:
%       ClearCC: connected component structure without CC on border

A = and(Frame,BW);
AIndx = find(A);
Indx = ones(CC.NumObjects,1);
for i=1:CC.NumObjects
    if any(ismember(AIndx,CC.PixelIdxList{i}))
        Indx(i) = 0;
    end
end
Indx = logical(Indx);

CC.PixelIdxList = CC.PixelIdxList(Indx);
CC.NumObjects = sum(Indx);
ClearCC = CC;

end

