
function [post, llt] = em_tc1_post(x, pp, mu, sig, pattern, xtb)
% compute the posterior z|x and the log-likelihood
%

N = size(x,1);
K = size(mu,1);

[post, ll] = em_censor_post(x, pp, mu, sig, pattern);

llt = ll - N*log(pp*(1-normcdf(xtb(1)*ones(K,1), mu(:,1), sqrt(squeeze(sig(1,1,:))))));

end
