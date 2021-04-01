function l = linePts(x1,x2,y1,y2)
%LINEPTS define line points between (x1,y1),(x2,y2)
% Input:
%       cordinate of two points involved in line
% output:
%       line points for all discrete points between two points

if abs(y2-y1) > abs(x2-x1)
    f = polyfit([y1,y2],[x1,x2],1);
    l = polyval(f,[y1:y2,y2:y1]);
else
    f = polyfit([x1,x2],[y1,y2],1);
    l = polyval(f,[x1:x2,x2:x1]);
end


end

