function [Idx,Samp] = circularSampling(O,R,Data,Eps)
%CIRCULARSAMPLING sampling data in a circle.
% Input:
%       O: [ 2 element vector (double)] center of circle
%       R: [scalar (double)] radius of circle
%       Data: [M x 2 (double)] data points in cartesian plane
%       Epx: [scalar(double)] Real circle radious = R + Epx (default = 0)
% Output:
%       Idx: [Z x 2 (logical)] index of qualified data
%       Samp: [Z x 2 (double)] data points inside of circle

if ~exist('Eps','var')
    Eps = 0;
end

% real circle radius
RR = R + Eps;

if isrow(Data)
    Data = Data';
end

Temp = Data;
Temp(:,1) = Temp(:,1) - O(1);
Temp(:,2) = Temp(:,2) - O(2);

Temp = sum(  Temp.^2  , 2 );

Idx = Temp <= RR.^2;

Samp = Data(Idx,:);

end

