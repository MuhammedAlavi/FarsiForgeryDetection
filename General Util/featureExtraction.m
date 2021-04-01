function FeatVec = featureExtraction(Tag, AllVar, ClustID, Weight)
%FEATUREEXTRACTION , reduce comparison matrix to a vector.
% Input:
%       Tag:[1 x 1 (logical)] Tag of comparison
%       AllVar: [N x M  (double)] Variables
%       ClustID: [vector (int)] id of cluster for Variable
%       Weight:[vector (double)] weight matrix
% Output:
%       FeatVec:[vector (double)] feature vector

% number of clusters
nClust = numel(unique(ClustID));

if isempty(Weight)
    Weight = ones(nClust,1);
end

% 3 features 1. mean 2. median 3.InterQuartileRange(IQR)
nFeature = 3;

% initialize feature vector
nVar = size(AllVar,2);
P = zeros(nClust, nVar * nFeature);


for i=1:nClust
    c = 1;
    I = ClustID == i;
    ClustVar = AllVar(I,:);
    for j=1:nVar
        P(i,c) = mean(ClustVar(:,j));
        c = c + 1;
        P(i,c) = median(ClustVar(:,j));
        c = c + 1;
        P(i,c) = iqr(ClustVar(:,j));
        c = c + 1;
    end
end

% weighted average for features
for i=1:nClust
    P(i,:) = P(i,:) .* Weight(i);
end
P = sum(P,1);
P = P ./ sum(Weight);

% concatenate tag at the end of feature vector
FeatVec = [P,Tag];
end

