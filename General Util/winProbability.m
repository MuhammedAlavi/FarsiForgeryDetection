function [WP] = winProbability( CmpMat,Ranking )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~ismatrix(CmpMat) && size(CmpMat,1) ~= size(CmpMat,2)
    error('input matrix is not square!');
end
N = size(CmpMat,1);

for i=1:N
    for j=i:N
        
        Iter = CmpMat(i,:);
    end
end
% remove
end

