function b = bernstein(n,t)
% BERNSTEIN: Bernstein Matrix
% b = bernstein(n,t)
% values b(j,k) of b^n_{k-1} at t(j)

% change to column vector (if necessary)
t=t(:);

% initialization for degree 0 
b = zeros(length(t),2*n+1); 
b(:,n+1) = 1;

% recursion for Bernstein polynomials
for k=1:n
   tk = repmat(t,1,2*n+1-k);
   b = (1-tk).*b(:,2:end) + tk.*b(:,1:end-1);
end
