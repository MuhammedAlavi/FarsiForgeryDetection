function BetterGray = grayEnhance(Gray,Low,High)
%GRAYENHANCE enhance gray level
% Input:
%       gray: gray level image
% Output:
%       betterGray: improved gray level

if ~exist('Low','var') || isempty(Low)
    Low = [.75,1];
end
if ~exist('High','var') || isempty(High)
    High = [0,1];
end
h = fspecial('average');
Gray = imfilter(Gray,h);
Gray = imadjust(Gray,Low,High);
BetterGray = Gray;

% crate a fake border, this line remove it
Plate = ones(size(BetterGray));
Plate = cast(Plate,class(BetterGray));
Plate = Plate .* intmax(class(BetterGray));
Plate(2:end-1,2:end-1) = 0;
BetterGray = BetterGray + Plate;
end

