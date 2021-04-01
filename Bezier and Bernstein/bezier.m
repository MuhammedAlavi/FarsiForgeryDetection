function Bez = bezier( X,Y,T,BezHelper,IsViz,varargin )
% BEZIER return bezier curve for control 
% points as input
%   Inputs: 
%       [x,y]: control points 
%       t: time point (default linspace(0,1,lengthOfCurvePoints))
%       IsViz: visualize the result (default = false)
%   Output: 
%       bez: bezier curve
%   ARBITRARY PARAMETERS:
%   arbitrary parameters are in this order
%   1. control point color
%   2. control point trajectory color
%   3. Bezier curve color
%   4. Bezier curve length
%   add these paramters at the end of function parameters
%   list with with just above order.
X = X(:);
Y = Y(:);
points = [X,Y];
all_pts = size(points,1);
assert(all_pts > 0,'empty points');

if ~exist('T','var') || isempty(T)
    T = linspace(0,1,all_pts);
end
if ~exist('IsViz','var')
    IsViz = false;
end
if exist('BezHelper','var') && ~isempty(BezHelper) && all_pts <= length(BezHelper)
    berMat = BezHelper{all_pts};
else
    berMat = bernstein(all_pts-1,T);
end
Bez = berMat * points;

%show the result?
params = {'black',30,2,'blue',2};%default values
for i=1:nargin-5  % for input parameters except 4 first one
    params{i} = varargin{i};
end
% initiated parameters
cntr_colr = params{1};
cntLiWid = params{3};
MK = params{2};
bez_colr = params{4};
LW = params{5};

if IsViz
    figure;
    hold on;
    scatter(X,Y,MK,cntr_colr,'fill');
    plot(X,Y,cntr_colr,'LineWidth',cntLiWid);
    scatter(Bez(:,1),Bez(:,2),MK,bez_colr,'fill');
    plot(Bez(:,1),Bez(:,2),bez_colr,'LineWidth',LW);
    % add the legand
    legend('Control Points','Main Trajectory',...
           'Bezier Points','Bezier Curve');
end
end

