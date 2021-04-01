function bez = SBezier(x,y,sliceTime,sliceLen,sliceOverlap,isViz,varargin)
%SIMPLEBEZIER generate a new bezier curve
% with breaking original curve to pieces and
% calculate bezier curve

% validate input arguments
validateattributes(sliceLen,{'numeric'},...
    {'odd','scalar','nonempty','finite'},mfilename,'sliceLen',4);
validateattributes(sliceOverlap,{'numeric'},...
    {'scalar','nonempty','finite','>=',0,'<',sliceLen},...
    mfilename,'sliceOverLap',5);
validateattributes(x,{'numeric'},...
    {'vector','nonempty','finite'},mfilename,'x',1);
validateattributes(y,{'numeric'},...
    {'vector','nonempty','finite'},mfilename,'y',2);

assert(size(x,1) == size(y,1));
curveLen = size(x,1);

% validate time of each slice
if ~exist('sliceTime','var') || isempty(sliceTime)
    sliceTime = linspace(0,1,sliceLen);
else
    validateattributes(sliceTime,{'numeric'},...
        {'vector','finite','numel',sliceLen,'>=',0,'<=',1},mfilename,'t',3);
end

% if not define the arguments
if ~exist('isViz','var')
    isViz = false;
end

% if length of curve less than or equal with slice
if curveLen <= sliceLen
    bez = bezier(x,y,sliceTime,isViz,varargin);
    return
end

% bernstein matrix
berMat = bernstein(sliceLen - 1,sliceTime);

% initialize curve slice matrix
% first and second column for sum of points (x,y)
% and third column for number of calculation
sliceMat = zeros(curveLen,3);

% how much is step for each iteration
step = sliceLen - sliceOverlap;
% calculate bezier of each segment
for i=1:step:curveLen
    if curveLen-sliceLen >= i
        % if segment is not at the
        % last part of curve
        try
        startpt = i;
        endpt = i + sliceLen - 1;
        cut = sliceMat(startpt:endpt,:);
        cut(:,1) = cut(:,1) + (berMat * x(startpt:endpt));
        cut(:,2) = cut(:,2) + (berMat * y(startpt:endpt));
        cut(:,3) = cut(:,3) + 1;
        catch
            disp('exc')
        end
    else
        % to calculate last part of curve
        startpt = curveLen - sliceLen + 1;
        endpt = curveLen;
        cut = sliceMat(startpt:endpt,:);
        cut(:,1) = cut(:,1) + (berMat * x(startpt:endpt));
        cut(:,2) = cut(:,2) + (berMat * y(startpt:endpt));
        cut(:,3) = cut(:,3) + 1;
    end
    sliceMat(startpt:endpt,:) = cut;
end

% calculate the average point
bez = zeros(curveLen,2);
bez(:,1) = sliceMat(:,1)./sliceMat(:,3);
bez(:,2) = sliceMat(:,2)./sliceMat(:,3);

%show the result?
params = {'black',30,2,'blue',2};%default values
for i=1:length(varargin)  % for input parameters except 4 first one
    params{i} = varargin{i};
end
% initiated parameters
cntr_colr = params{1};
cntLiWid = params{3};
MK = params{2};
bez_colr = params{4};
LW = params{5};

if isViz
    figure;
    hold on;
    scatter(x,y,MK,cntr_colr,'fill');
    plot(x,y,cntr_colr,'LineWidth',cntLiWid);
    scatter(bez(:,1),bez(:,2),MK,bez_colr,'fill');
    plot(bez(:,1),bez(:,2),bez_colr,'LineWidth',LW);
    % add the legand
    legend('Control Points','Main Trajectory',...
           'Bezier Points','Bezier Curve');
end
end