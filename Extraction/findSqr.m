function [CC,SqrImg] = findSqr(HandWritImg,SqrSz)
%FINDSQR search and find four squares that indicates 
% four angle of handwriting
% Input:
%       HandWritImg: black and white image (black is true)
%       SqrSz: size of Square
% Output:
%       CC: connected component structure of image
%       SqrImg: image with 4 square indicator

ImgCC = bwconncomp(HandWritImg);
SqrElements = SqrSz^2;
Index = false(ImgCC.NumObjects,1);

% search all cc
for i=1:ImgCC.NumObjects
    if numel(ImgCC.PixelIdxList{i}) == SqrElements
        % convert index CC to 2D image
        ImageOfCC = ind2img(ImgCC.PixelIdxList{i},HandWritImg);
        
        % if selected connected component is exactly a square
        if size(ImageOfCC,1) == size(ImageOfCC,2)
            Index(i) = true;
        end
    end
end

% select connected components
CC = ImgCC;
CC.PixelIdxList = CC.PixelIdxList(Index);
CC.NumObjects = sum(Index);

% if user needs image of 4 squares
if nargout > 1
    SqrImg = false(size(HandWritImg));
    for i=1:CC.NumObjects
        SqrImg(CC.PixelIdxList{i}) = true;
    end
end

end

