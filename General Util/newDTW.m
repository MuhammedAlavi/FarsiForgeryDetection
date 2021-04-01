function D = newDTW(X,Y,WindowSz,XLen,YLen)
% dtw distance

if ~exist('TargetLen','var')    
    XLen = length(X);
end

if ~exist('TestLen','var')
    YLen = length(Y);
end

DTW = inf(XLen+1,YLen+1);

% preallocation of measurement matrix
DTW(1,1) = 0;

% define window size for dtw
D = XLen-YLen;

% absolute value of D
if D<0
    D = -D;
end

% window size works if => abs(TargetLen,TestLen) < window size
if ~exist('WindowSz','var')
    WindowSz = YLen;
elseif isempty(WindowSz)
    WindowSz = YLen;
elseif D > WindowSz
    WindowSz = D;
end

% create measurement matrix
for i=2:XLen+1
    % lower bound
    LB = 2;
    if i-WindowSz > LB
        LB = i-WindowSz;
    end
    
    % upper bound
    UB = YLen+1;
    if i+WindowSz < UB
        UB = i+WindowSz;
    end
    
    for j=LB:UB
        Cost = X(i-1) - Y(j-1);
        
        % absolute value of Cost
        if Cost < 0
            Cost = -Cost;
        end
        
        % check three adjacent elements of measurement matrix
        A = DTW(i-1,j-1);
        B = DTW(i-1,j);
        C = DTW(i,j-1);
        
        if A <= B && A <= C
            DTW(i,j) = Cost + A;
        elseif B < A && B < C
            DTW(i,j) = Cost + B;
        elseif C < A && C < B
            DTW(i,j) = Cost + C;
        end
    end
end
D = DTW(XLen+1,YLen+1);


end

