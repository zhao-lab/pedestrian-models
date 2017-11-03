function param = em_tc1(x, K, xtb, xcb, init)
% gaussian mixture model em algorithm for censored data
% 
% x         N*d     censored observations
% K                 number of components
% xtb       2*d     truncation bounds [lower ; upper]
% xcb       2*d     censoring bounds [lower ; upper]
% init      struct  initial parameters
% 

if nargin<3 || ~exist('xtb', 'var')
    xtb = [-Inf; Inf];
end

if nargin<4 || ~exist('xcb', 'var')
    xcb = [-Inf; Inf];
end


% const
[N,d] = size(x);

regVal = eps(max(eig(cov(x))));        % a small number added to diagonal of cov mtx
tol = 1e-10;            % termination tolerance
maxIt = 100*d;          % max number of iterations


% cencoring pattern
pattern = em_censor_pattern(x, xcb);


% init
if nargin < 5,  init = struct([]);  end
param0 = em_init(x, K, init, pattern);
pp = param0.pp;
mu = param0.mu;
sig = param0.C;


% 
ll_old = -Inf;
ll_hist = zeros(maxIt,1);

% EM
for it = 1:maxIt
    
    % E
    % compute the posterior z|x and log-likelihood
    [post, ll] = em_tc1_post(x, pp, mu, sig, pattern, xtb);
    
    
    % check convergence
    lldiff = ll - ll_old;
    ll_hist(it) = ll;
    if (lldiff >= 0 && lldiff < tol*abs(ll))
        break;
    end
    ll_old = ll;
    
    
    % unnormalized truncated component probability
    ppd = sum(post);
    
    
    % M
    for k=1:K
        % compute sufficient statistics
        [xhat_e, Q_e, alpha] = xhatQ(x, mu(k,:), sig(:,:,k), pattern);
        
        % Numerical issue
        % when alpha is below the machine precision, 
        % alpha=0, xhat=inf and Q=nan
        alpha0 = (alpha==0);
        xhat = xhat_e;
        xhat(alpha0,:) = 0;
        Q = Q_e;
        Q(:,:,alpha0) = 0;
        
        
        % mu
        % correct the truncation bias
        % TODO : truncation is only on the 1st coordinate
        mk = normpdf(xtb(1), mu(k,1), sqrt(sig(1,1,k))) ...
           / (1 - normcdf(xtb(1), mu(k,1), sqrt(sig(1,1,k))));
        
        mu(k,:) = (post(:,k)'*xhat) / ppd(k) ...
            - sig(1,:,k) * mk;
        
        
        % sig
        % correct the truncation bias
        % TODO : truncation is only on the 1st coordinate
        Rk = (xtb(1) - mu(k,1)) * normpdf(xtb(1), mu(k,1), sqrt(sig(1,1,k))) ...
            / sig(1,1,k) ...
            / (1 - normcdf(xtb(1), mu(k,1), sqrt(sig(1,1,k))));
        
        xhat0 = bsxfun(@minus, xhat, mu(k,:));
        xhat0 = bsxfun(@times, sqrt(post(:,k)), xhat0);
        
        signew = xhat0'*xhat0 + reshape(reshape(Q, [d*d,N]) * post(:,k), [d,d]);
        signew = signew / ppd(k) ...
            - sig(1,:,k)'*sig(1,:,k)*Rk;
        
        sig(:,:,k) = (signew+signew')/2 + regVal*eye(d);
    end
    
    
    % recover component probability
    % TODO : truncation is only on the 1st coordinate
    ppdratio = 1 - normcdf(xtb(1)*ones(K,1), mu(:,1), sqrt(squeeze(sig(1,1,:))));
    pp = ppd ./ ppdratio';
    
    
    % normalize component probability
    pp = pp / sum(pp);
    
end


% # of parameters
nParam = (d-1) + K*(d + d*(d+1)/2);


% output
param.K = K;
param.pp = pp;
param.mu = mu;
param.C = sig;

param.iters = it;
param.log_lh = ll;
param.AIC = -2*ll + 2*nParam;
param.BIC = -2*ll + nParam*log(N);
param.ll_hist = ll_hist(1:it);
param.regVal = regVal;

end



function [xhat, Q, alpha] = xhatQ(x, mu, sig, pattern)

[N,d] = size(x);
x0 = bsxfun(@minus, x, mu);

xhat = x;
Q = zeros(d,d,N);
alpha = ones(N,1);

pttn = pattern.uniq_censored;
cnt = pattern.count;
xpttn = pattern.xpttn;
xR = pattern.range;

for ii=1:length(cnt)
    on = pttn(ii,:);
    mn = ~on;
%     dm = sum(mn);
    idx = find(xpttn(:,ii));
    
    % check censored element
    if all(on)
        continue
    end
    
    [R,f] = chol(sig(on,on));
    
    % mean and covariance of conditional normal 
    mu_mo = (x0(idx,on) / R) * (sig(mn,on) / R)';
    mu_mo = bsxfun(@plus, mu_mo, mu(mn));
    
    sigRinv = sig(mn,on) / R;
    sig_mo = sig(mn,mn) - sigRinv*sigRinv';
    
    % mean and covariance of truncated normal 
    xRidx = xR(idx,:);
    xRl = cell2mat(xRidx(:,1));
    xRu = cell2mat(xRidx(:,2));
    [tmu, tcov, talpha] = tmvn_m3(mu_mo, sig_mo, xRl, xRu);
    
    
    % replace censored elements with conditional mean
    xhat(idx,mn) = tmu;
    % covariance corrections for censored elements
    Q(mn,mn,idx) = tcov;
    % probability mass in the truncated region
    alpha(idx) = talpha;
end

end


function param0 = em_init(x, K, init, pattern)

[~,d] = size(x);

cmpl = pattern.complete;

if isempty(init)
    xc = x(cmpl,:);
    [idx, cent] = kmeans(xc, K);
    
    pp = zeros(1,K);
    mu = zeros(K,d);
    sig= zeros(d,d,K);
    for k=1:K
        pp(k) = sum(idx==k);
        mu(k,:) = cent(k,:);
        if pp(k)
            sig(:,:,k) = cov(xc(idx==k,:));
        else
            sig(:,:,k) = eye(d)*1e-10;
        end
    end
    
    param0.pp = pp;
    param0.mu = mu;
    param0.C = sig;
else
    param0 = init;
end


end

