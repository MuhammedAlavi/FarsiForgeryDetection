function CCImg = ind2img(Index,Img)
%IND2IMG converts index connected component to image
% Input:
%       Index: index of connected component
%       ImgSz: size of image
% Output:
%       Img: image of connected component

% check type of image
if isinteger(Img)
    c = class(Img);        
    if ndims(Img) == 3
        % create a blank white page
        CCImg = ones(size(Img),c);
        CCImg = CCImg.* intmax(c);
        % --------------------
        R =  Img(:,:,1);
        G =  Img(:,:,2);
        B =  Img(:,:,3);
        % --------------------
        Temp_R = CCImg(:,:,1);
        Temp_G = CCImg(:,:,2);
        Temp_B = CCImg(:,:,3);
        % --------------------
        Temp_R(Index) = R(Index);
        Temp_G(Index) = G(Index);
        Temp_B(Index) = B(Index);
        % --------------------
        CCImg = cat(3,Temp_R,Temp_G,Temp_B);
        
    elseif ismatrix(Img)
        % if image is a gray scale
        CCImg = zeros(size(Img),c);
        CCImg(Index) = Img(Index);
    end

end

if islogical(Img)
    CCImg = false(size(Img));
    CCImg(Index) = 1;
end

[Row,Col] = ind2sub([size(Img,1),size(Img,2)],Index);
s = [min(Row),min(Col)];
width = max(Col)-min(Col);
len = max(Row)-min(Row);

%crop image
if isinteger(CCImg)
    if ndims(CCImg) == 3
        CCImg = CCImg(s(1):s(1)+len,s(2):s(2)+width,1:3);
    elseif ismatrix(CCImg)
        CCImg = CCImg(s(1):s(1)+len,s(2):s(2)+width);
    end
elseif islogical(CCImg)
    CCImg = CCImg(s(1):s(1)+len,s(2):s(2)+width);
end
end
