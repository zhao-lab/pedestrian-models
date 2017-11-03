function [x, z] = util_gmm(pi, mean, sigma, N)
% generate gmm dataset
% 
%   pi      1 x K       K component weights
%   mean    d x K
%   sigma   d x d x K   covariance matrix
% 
%   x       N x d
%   z       N x 1
% 
if (size(pi,1)>1), pi = pi'; end
if abs(sum(pi)-1)>eps, error('pi is not proper'); end

% component
K = length(pi);

% component weight
pi = [0 cumsum(pi(1:end-1))];

% latent variable
z = sum(repmat(rand(N,1), 1, K) >= repmat(pi, N, 1), 2);

% dataset
x = zeros(N,size(mean,1));
for k=1:K
    x(z==k,:) = mvnrnd(mean(:,k)', sigma(:,:,k), sum(z==k));
end

if K<=2, z = logical(z-1); end
