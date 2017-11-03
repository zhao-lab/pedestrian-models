
function [post, ll, logpdf] = em_std_post(x, pp, mu, sig)

N = size(x,1);
K = size(mu,1);

log_lh = zeros(N, K);

for k=1:K
    log_lh(:,k) = util_lnpost(x, pp(k), mu(k,:), sig(:,:,k));
end

[post, ll, logpdf] = util_post(log_lh);

end


function log_lh = util_lnpost(x, p, mu, sigma)

log_p = log(p);
[n,d] = size(x);
x0 = bsxfun(@minus, x, mu);
[R,f] = chol(sigma);
xRinv = x0 / R;
quadform = sum(xRinv.^2, 2);
logSqrtDetSig = sum(log(diag(R)));
log_lh = -0.5*quadform + (-logSqrtDetSig + log_p) - d*log(2*pi)/2;

end % function util_lnpost


function [post, ll, logpdf] = util_post(log_lh)

maxll = max(log_lh,[],2);
post = exp(bsxfun(@minus, log_lh, maxll));
density = sum(post,2);
post = bsxfun(@rdivide, post, density);
logpdf = log(density) + maxll;
ll = sum(logpdf);

end % function util_post

