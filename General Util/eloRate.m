function ExpectedScore = eloRate(Ra,Rb,K)
%ELORATE calculate expected score according to elo rating.
if ~exist('K','var')
    K = 400;
end
Qa = 10^(Ra/K);
Qb = 10^(Rb/K);

ExpectedScore = Qa/(Qa+Qb);

end

