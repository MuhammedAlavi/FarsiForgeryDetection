function [r,C] = curvature(P,Q,R)
%CURVATURE is a circle passed on P,Q,R. P,Q are two side points and Q is
% middle point.
% Input:
%       P,Q,R: [vector with 2 elements]
%       first element is X and second is Y
% Output:
%       r: [scalar] radious of circle

P = P(:);
Q = Q(:);
R = R(:);

% first chord
d1 = Q-P;
mid1 = (P + Q) ./ 2;
x01 = mid1(1);
y01 = mid1(2);
m1 = d1(2)/d1(1);
m1 = -1/m1;

% second chord
d2 = R-Q;
mid2 = (R + Q) ./ 2 ;
x02 = mid2(1);
y02 = mid2(2);
m2 = d2(2)/d2(1);
m2 = -1/m2;

% common point
x = (y01 - (m1 * x01) - y02 + (m2 * x02)) / (m2 - m1);
y = m1 * (x - x01) + y01;

% % % % figure,hold on
% % % % scatter(P(1),P(2),'fill')
% % % % scatter(Q(1),Q(2),'fill')
% % % % scatter(R(1),R(2),'fill')
% % % % scatter(x01,y01,'fill')
% % % % scatter(x02,y02,'fill')
% % % % scatter(x,y,'fill','black')
% % % % plot([P(1),Q(1)],[P(2),Q(2)]);
% % % % plot([R(1),Q(1)],[R(2),Q(2)]);
% % % % plot([mid1(1),x],[mid1(2),y])
% % % % plot([mid2(1),x],[mid2(2),y])

% % % v1 = mid2 - R;
% % % v2 = mid2 - [x;y];
% % % ab = sqrt(sum(v1.^2))*sqrt(sum(v2.^2));
% % % theta = acosd(sum(v1.*v2)/ab);



% radious of circle
C = [x,y];
D = Q - C';
r = sqrt(sum(D.^2));


end

