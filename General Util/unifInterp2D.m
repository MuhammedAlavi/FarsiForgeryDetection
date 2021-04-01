function P = unifInterp2D(P,K)
%UNIFINTERP2D interpolate uniformly with divide each consequent points by
% two K times.
for i=1:K
    v = diff(P);
    v = v ./ 2;
    v = v + P(1:end-1,:);
    NewP = zeros((size(P,1)+size(v,1)),2);
    NewP(1:2:end,:) = P;
    NewP(2:2:end,:) = v;
    P = NewP;
end

end

