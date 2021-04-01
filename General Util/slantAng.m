function Ang = slantAng(Img)
%SLANTANG angle of slant of CC
% Input:
%       Img:[N x M (logical)] black and white image
% Output:
%       Ang:[scalar(double)] angle of slant(degree)
% Algorithm from:
% Cheriet, Mohamed, et al. Character recognition systems:
% a guide for students and practitioners.
% John Wiley & Sons, 2007.

[R,C] = find(Img);
[~,i] = min(R+C);
Beta2 = sum([R(i),C(i)]);

[~,i] = max(R+C);
Beta1 = sum([R(i),C(i)]);


[~,i] = min(R - C);
Beta4 = sum([R(i),C(i)]);

[~,i] = max(R - C);
Beta3 = sum([R(i),C(i)]);


% rows of image
I = sum(Img,2);
I = I > 0;
A = abs(find(I,1,'first')-find(I,1,'last'));

B = (Beta4 + Beta1 - Beta3-Beta2) / 2;

Ang = atan2d(A,B);
end

