function NewPts = slantCorrection(X,Y,Center,Ang)
%SLANTCORRECTION, transfer points to (0,0) and rotate bezier
% curve around center of mass.
% Input:
%       Ang: [1 x 1 (int)] degree of rotation
%       +: rotate counterclockwise
%       -: colckwise
% Output:
%       new points transfered to (0,0) and rotated

X = X(:);
Y = Y(:);

% transfer center of mass to (0,0)
X = X - Center(1);
Y = Y - Center(2);

C = cosd(Ang);
S = sind(Ang);

% rotation matrix
% |cos(theta)  -sin(theta)|
% |sin(theta)   cos(theta)|
NewX = X.*C - Y.*S;
NewY = X.*S + Y.*C;

% return back points
NewX = NewX + Center(1);
NewY = NewY + Center(2);

NewPts = [NewX,NewY];


end

