
function [post, ll, logl] = em_censor_post(x, pp, mu, sig, pattern)
% compute the posterior z|x and the log-likelihood
% 

N = size(x,1);
K = size(mu,1);

log_lh = zeros(N,K);

for k=1:K
    log_lh(:,k) = em_lnpost(x, pp(k), mu(k,:), sig(:,:,k), pattern);
end

[post, ll, logl] = em_post(log_lh);

end



function log_lh = em_lnpost(x, pp, mu, sig, pattern)
% unnormalized weighted log likelihood of x
% 

N = size(x,1);

log_p = log(pp);    % weight : log of component probability
x0 = bsxfun(@minus, x, mu);

log_lh = zeros(N,1);

pttn = pattern.uniq_censored;
cnt = pattern.count;
xpttn = pattern.xpttn;
xR = pattern.range;

for ii=1:length(cnt)
    on = pttn(ii,:);
    mn = ~on;
    do = sum(on);
    dm = sum(mn);
    idx = find(xpttn(:,ii));
    
    % log likelihhod of observed x_o
    [R,f] = chol(sig(on,on));
    xRinv = x0(idx, on) / R;
    quadform = sum(xRinv.^2, 2);
    logSqrtDetSig = sum(log(diag(R)));
    
    % log prob mass of censored x_m given x_o
    if any(mn)
        mu_mo = (x0(idx,on) / R) * (sig(mn,on) / R)';
        mu_mo = bsxfun(@plus, mu_mo, mu(mn));
        
        sigRinv = sig(mn,on) / R;
        sig_mo = sig(mn,mn) - sigRinv*sigRinv';
        
        %%%
        varless0 = diag(sig_mo) < eps;
        if any(varless0)
            sig_mo(diag(varless0)) = eps;
        end
        
        xRidx = xR(idx,:);
        xRl = cell2mat(xRidx(:,1));
        xRu = cell2mat(xRidx(:,2));
        log_Phi = log(mvncdf(xRl - mu_mo, xRu - mu_mo, zeros(1,dm), sig_mo));
    else
        log_Phi = 0;
    end
    
    % weighted log likelihood of x
    log_lh(idx) = -0.5*quadform + (-logSqrtDetSig + log_p + log_Phi) - do*log(2*pi)/2;
end

end



function [post, ll, logl] = em_post(log_lh)
% normalize the weighted log likelihood of x
% 

maxll = max(log_lh,[],2);
post = exp(bsxfun(@minus, log_lh, maxll));
density = sum(post,2);
post = bsxfun(@rdivide, post, density);
logl = log(density) + maxll;
ll = sum(logl);

end

