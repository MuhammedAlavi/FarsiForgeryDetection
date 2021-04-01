function Fun = normalInterp(Y,N,Method)
%NORMALINTERP, normalize points Xq between [0,1] in N points.
% Input:
%       Y: [vector (double)] points of curve
%       N: [1 x 1 (int)] number of points in interval [0,1] for interpolation
% Output:
%       Yq: [vector (double)] interpolated points

if ~exist('Method','var')
    Method = 'linear';
end
Y = Y(:)';
Len = length(Y);
Xv = linspace(0,1,Len);
Xv = Xv(:)';
Xq = linspace(0,1,N);
Xq = Xq(:)';
Yq = interp1(Xv,Y,Xq,Method);
Xq = Xq(:);
Yq = Yq(:);
Fun = [Xq,Yq];
end

