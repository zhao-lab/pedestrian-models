function param = em_std(x, K, init)
global aic bic
str_start.mu = init.mu;
str_start.Sigma = init.C;
str_start.PComponents = init.pp;
options = statset('MaxIter',2000);
%obj = gmdistribution.fit(X,2,'Options',options);

%res_obj = fitgmdist(x, K, 'Start', str_start,'Options',options);
res_obj = fitgmdist(x, K, 'Options',options);
aic(K) = res_obj.AIC;
bic(K) = res_obj.BIC;
param.pp = res_obj.PComponents;
param.mu = res_obj.mu;
param.C = res_obj.Sigma;
