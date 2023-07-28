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

% fit the model 
b  = (pinv(X*iC*X')*X*iC)'*y;

% iterative fit
% func = @(b) X*spm_unvec(b,b0);
% Y  = y;
% x0 = zeros(size(X,2),size(y,2)) + 1./length(x0(:));
% x0 = x0(:);
% V  = x0 + 1/8;
% 
% 
% op = AO('options');
% op.fun = func;       % function/model f(x0)
% op.x0  = x0(:);      % start values: x0
% op.y   = Y(:);       % data we're fitting (for computation of objective fun, e.g. e = Y - f(x)
% op.V   = V(:);       % variance / step for each parameter, e.g. ones(length(x0),1)/8
% op.hyperparams=0;
% op.objective='gauss'; % select smooth Gaussian error function
% 
% [X,F,CV,~,Hi] = AO(op); 


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