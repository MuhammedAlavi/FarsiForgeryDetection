function [dist,ix,iy] = DTW(x,y,metric) 
%DTW    Distance between signals via Dynamic Time Warping
%   DIST = DTW(X,Y) computes the minimum Euclidean distance, DIST, between
%   two vectors X and Y, where samples in either X or Y may consecutively
%   repeat any number of times.  If X and Y are matrices, then X and Y must
%   have the same number of rows and DTW will allow columns of X and Y to
%   consecutively repeat.
%
%   [DIST,IX,IY] = DTW(X,Y) additionally returns the warping path, IX and
%   IY, that minimizes the total Euclidean distance between X(IX) and Y(IY)
%   when X and Y are vectors and between X(:,IX) and Y(:,IY) when X and Y
%   are matrices.
%
%   [DIST,IX,IY] = DTW(X,Y,MAXSAMP) additionally restricts IX and IY so
%   that they must be within MAXSAMP samples of a straight-line fit between
%   X and Y.  If MAXSAMP is unspecified, then no restriction will be placed
%   upon IX or IY.
%
%   [DIST,IX,IY] = DTW(X,Y,...,METRIC) will return in DIST the summed
%   distance between the corresponding entries of X and Y according to the
%   distance metric defined by METRIC.  The default is 'euclidean'.
%      'absolute'  the sum of the absolute (manhattan) differences
%      'euclidean' the root sum squared differences
%      'squared'   the squared Euclidean distance
%      'symmkl'    the symmetric Kullback-Leibler distance 
%
%   DTW(...) without output arguments plots the original and aligned
%   signals.  DTW displays the alignment of X and Y via a line plot when
%   they are vectors, and via horizontally aligned images when they are
%   matrices.  If the matrices are complex, the real and imaginary
%   portions appear in the top and bottom half of each image, respectively.
%   
%   % Example 1:
%   %   Compute and plot the best Euclidean distance between real
%   %   chirp and sinusoidal signals using dynamic time warping.
%   x = chirp(0:999,0,1000,1/100);
%   y = cos(2*pi*5*(0:199)/200);
%   dtw(x,y)
%
%   % Example 2:
%   %   Compute and plot the best Euclidean distance between complex
%   %   chirp and sinusoidal signals using dynamic time warping.
%   x = exp(2i*pi*(3*(1:1000)/1000).^2);
%   y = exp(2i*pi*9*(1:399)/400);
%   dtw(x,y)
%
%   % Example 3:
%   %   Align handwriting samples along the x-axis.
%   x = double(imread('MATLAB1.gif'));
%   y = double(imread('MATLAB2.gif'));
%   dtw(x,y);
%
%   See also ALIGNSIGNALS, FINDDELAY, XCORR.

%   References: 
%   * H. Sakoe and S. Chiba, "Dynamic Programming Algorithm Optimization
%     for Spoken Word Recognition" IEEE Transactions on Acoustics, Speech
%     and Signal Processing, Vol. ASSP-26, No. 1, Feb 1978, pp. 43-49.
%   * K.K. Paliwal, A. Agarwal and S.S. Sinha, "A Modification over Sakoe
%     and Chiba's dynamic time warping algorithm for isolated word
%     recognition" IEEE International Conference on ICASSP 1982., Vol. 7,
%     pp. 1259-1261

%   Copyright 2015 The MathWorks, Inc.

if ~exist('metric','var')
    metric = 'euclidean';
end

x = x(:)';
y = y(:)';
C = dtwmex(x, y, metric);
dist=C(size(x,2),size(y,2));
[ix,iy] = traceback(C);

end
%-------------------------------------------------------------------------
function [ix,iy] = traceback(C)
m = size(C,1);
n = size(C,2);

% pre-allocate to the maximum warping path size.
ix = zeros(m+n,1);
iy = zeros(m+n,1);

ix(1) = m;
iy(1) = n;

i = m;
j = n;
k = 1;

while i>1 || j>1
  if j == 1
    i = i-1;
  elseif i == 1
    j = j-1;
  else
    % trace back to the origin, ignoring any NaN value
    % prefer i in a tie between i and j
    cij = C(i-1,j-1);
    ci = C(i-1,j);
    cj = C(i,j-1);
    i = i - (ci<=cj | cij<=cj | cj~=cj);
    j = j - (cj<ci | cij<=ci | ci~=ci);
  end
  k = k+1;
  ix(k) = i;
  iy(k) = j;
end

ix = ix(k:-1:1);
iy = iy(k:-1:1);
end