function [Radius,Center] = curveCurvature(Curve,ChordLen)
%curveCurvature calculate curvature of a curve with different chord length.

Len = length(Curve);
First = 1;
Mid = 1 + ChordLen;
Last = 1 + 2*ChordLen;

Radius = zeros(Len,1);
Center = zeros(Len,2);
i = 0;
while Last <= Len
    P = Curve(First,:);
    Q = Curve(Mid,:);
    R = Curve(Last,:);
    i = i + 1;
    [Radius(i),Center(i,:)] = curvature(P,Q,R);
    First = First + 1;
    Mid = Mid + 1;
    Last = Last + 1;
end

% remove zero elements
Radius = Radius(Radius > 0);

end

