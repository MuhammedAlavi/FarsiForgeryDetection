function [DDirection] = connCompSearch(TargetCC,CCList,Z,W)
%CONNCOMPCMP return result of Z most resemble of TargetCC with CCList.
% Input:
%       TargetCC:[ConnComp object] one object of CC
%       CCList:[ConnComp object] list of object of CC to search
%       Z: [scalar (int)] number of best CC to choose
%       W: [scalar (int)] size of window for dtw algorithm
% Output:
%       CmpResult: [Matrix] comparison matrix
% compare result include:

if ~isa(TargetCC,'ConnComp')
    error('Target Object must be a ConnComp!');
end

if ~isa(CCList(1),'ConnComp')
    error('Search List must be a ConnComp object!');
end

if isempty(CCList) || isempty(TargetCC)
    DDirection = [];
    return
end
if ~exist('Z','var')
    Z = 1;
end

% number of CCList elements
N = numel(CCList);

DDirection = zeros(N,11);

TargetAng = TargetCC.Ang;
TargetID = TargetCC.ID;
TargetInkModel = TargetCC.IM;
TargetBProjHist = TargetCC.BProjHist;
TargetSProjHist = TargetCC.SProjHist;
TargetCOMProp= TargetCC.COMProp;
TargetSlant = TargetCC.Slant;

for i=1:N
  % calculate dtw of ink model and direction vector
    DDirection(i,1) = TargetID;
    DDirection(i,2) = CCList(i).ID;
    DDirection(i,3) = dtw_c(TargetAng,CCList(i).Ang,W);
    DDirection(i,3) = DDirection(i,3) / length(CCList(i).Ang);
end

DD = DDirection(:,3);

% sort distance and return Z closest points
[~,I] = sort(DD,'ascend');

if Z < numel(I)
    I = I(1:Z);   
end

% choose best z points
DDirection = DDirection(I,:);

if Z > numel(I)
    l = numel(I);
else
    l = Z;
end

for i=1:l
    IMLen = length(TargetInkModel);
    
    % check an exception conditon
    if length(TargetInkModel) == 2
        TargetInkWithAng = abs(TargetAng) ./ TargetInkModel(end) ;
    else
        TargetInkWithAng = abs(diff(TargetAng)) ./ TargetInkModel(2:end-1) ;
    end
    TestInkWithAng = CCList(I(i)).IM(1:end-1) ./ CCList(I(i)).Ang;
    
    DDirection(i,4) = dtw_c(TargetInkModel,CCList(I(i)).IM)/IMLen;
    DDirection(i,5) = dtw_c(TargetInkWithAng, TestInkWithAng);
    DDirection(i,6) = mean(TargetInkModel) - mean(CCList(I(i)).IM);
    DDirection(i,7) = dtw_c(TargetBProjHist,CCList(I(i)).BProjHist);
    DDirection(i,8) = dtw_c(TargetSProjHist,CCList(I(i)).SProjHist);
    DDirection(i,9) = TargetSlant - CCList(I(i)).Slant;
    DDirection(i,10:11) = TargetCOMProp - CCList(I(i)).COMProp;
%     ThumbSize = CCList(I(i)).ThumbSize;
%     ThumbSize = ThumbSize - TargetThumbSize;
%     ThumbSize = ThumbSize .* F;
%     DDirection(i,5) = ThumbSize(1);
%     DDirection(i,6) = ThumbSize(2);
end
    
% remove Target and Test CC identifiers 
DDirection = DDirection(:,3:end);
end

