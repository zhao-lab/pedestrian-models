function K = util_distM(x, y, mtx)
% compute pairwise Mahalanobis distance between x and y
% 
%   x       n x p feature vector
%   y       m x p feature vector
%   mtx     p x p covariance matrix
% 
%   K       n x m distant matrix
% 
if nargin < 2, y = x; end

[ry, cy] = size(y);
[rx, cx] = size(x);

if cx~=cy
    error('util_distM:dimension', ...
          'dim(x) ~= dim(y)');
end

if nargin < 3, mtx = eye(p); end


[R,f] = chol(mtx);
% invcov = inv(mtx);

K = zeros(rx, ry);
for ii=1:rx
    dd = repmat(x(ii,:),ry,1) - y;
%     ds = sum((dd*invcov).*dd,2);
    ddRinv = dd / R;
    ds = sum(ddRinv.^2, 2);
    K(ii,:) = ds';
end
