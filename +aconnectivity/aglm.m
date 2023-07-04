function [b,G,stats] = aglm(X,y)
% fit a multivariate GLM using the weighted least square estimator, where
% the weights are the pseudo inverse of the design matrix:
%
% [b,model,stats] = aglm(X,Y)
%
% implements: b  = (pinv(X*iC*X')*X*iC)'*y;
%
% returns betas, full model (X*b) and statistics including R2, F-value and
% p-value for full model.
%
% Student t-tests are used to assess each individual predictor, assumed
% within-subject (paired t).
%
% AS2023


% weighted least squares estimator - is the minimum-variance unbiased
% estimator (Gauss-Markov theorem)

% computes weights
iC = pinv(cov(X));iC(1,:) = 1; iC(:,1) = 1;

% fir the model ('\' picks appropriate algorithm?)
b  = (pinv(X*iC*X')*X*iC)'*y;

% full model
G = X*b;
r = y - G;

r2 = abs(1 - (norm(r)/norm(y-mean(y)))^2);

% convert to F
df1 = size(b,1) - 1;
df2 = abs(size(b,2) - (size(b,1)+1));

% stats on the GLM:
F = (r2 ./ (1-r2) ) .* (df2./df1);
F = abs(F);
p = 1-fcdf(F,df1,df2);

% stats to return
stats.e  = r;
stats.r2 = r2;
stats.F  = F;
stats.p  = p;

% compute stats (t) for each variable
[H,Xp,CI,STATS] = ttest(b');

stats.Xpval = Xp;
stats.Xtval = STATS.tstat;