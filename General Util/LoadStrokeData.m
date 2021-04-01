function D = LoadStrokeData(d)
%LOADSTROKEDATA, loads data for find a specific stroke.

switch d
    case 'Y'
        D = load('Y.mat');
        D = D.Y;
    case 'N'
        D = load('N.mat');
        D = D.N;
    otherwise
        D = [];
end
        


end

