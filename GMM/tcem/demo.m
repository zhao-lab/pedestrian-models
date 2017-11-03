% % % % % 
% demo

clear
home


% 
% data
N = 1000;

pp = [5/10, 2/10, 3/10];
mu1 = [-3, 3];
mu2 = [10, -1];
mu3 = [20, 20];
C1 = diag([20, 5]);
C2 = diag([5, 20]);
C3 = diag([20, 20]);
paramt.pp = pp;
paramt.mu = [mu1; mu2; mu3];
paramt.C = reshape([C1(:), C2(:), C3(:)], [2 2 3]);
clear pp mu1 mu2 mu3 C1 C2 C3

y = util_gmm(paramt.pp, paramt.mu', paramt.C, N);

% truncation and censoring
tc = 25;
idxt = y(:,1)>0;    % truncation
idxcu = y>tc;       % censoring
idxcl = y<0;        % censoring
idxc = any(idxcu | idxcl,2);

x = y;
x(idxcu) = tc;
x(idxcl) = 0;
x = x(idxt,:);
idxcx = idxc(idxt);

pattern = em_censor_pattern(x, [0; tc]);
xc = x(pattern.complete,:);


% 
% model fitting
K = 3;


% initialize by k-means 
%   run k-means five times and initialize each algorithm with the
%   parameters that achieves the highest likelihood
par_k = init_kmeans(x, K, 5);


% standard EM
par_k_std_ll = zeros(5,1);
for m=1:5
    [temp,par_k_std_ll(m)] = em_std_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C);
end
[temp,par_k_std_idx] = max(par_k_std_ll);
param1_init = par_k{par_k_std_idx};

param1 = em_std(x, K, param1_init);

[temp,idx] = sort(param1.mu(:,1), 'ascend');
param1.pp = param1.pp(idx);
param1.mu = param1.mu(idx,:);
param1.C = param1.C(:,:,idx);


% truncated and censored data EM
par_k_tc_ll = zeros(5,1);
for m=1:5
    [temp,par_k_tc_ll(m)] = em_tc1_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C, pattern, [0;Inf]);
end
[temp,par_k_tc_idx] = max(par_k_tc_ll);
param2_init = par_k{par_k_tc_idx};

param2 = em_tc1(x, K, [0;Inf], [0; tc], param2_init);

[temp,idx] = sort(param2.mu(:,1), 'ascend');
param2.pp = param2.pp(idx);
param2.mu = param2.mu(idx,:);
param2.C = param2.C(:,:,idx);


% 
% show the result
figure(2), clf
set(gcf, 'Position', [200 400 600 250])

ngrid = 30;
[xx1, xx2] = ndgrid(linspace(min(paramt.mu(:,1)),tc,ngrid));
xg = [xx1(:), xx2(:)];


% standard EM
subplot(121)
plot(x(:,1), x(:,2), 'xr', ...
    'markersize', 5)
hold on
plot(x(idxcx,1), x(idxcx,2), 'xr', ...
    'markersize', 5)
hold off

hold on
plot(paramt.mu(:,1), paramt.mu(:,2), '+k', ...
    param1.mu(:,1), param1.mu(:,2), 'ok', 'linewidth', 2)
hold off

for k=1:size(paramt.mu,1)
    xxd = util_distM(xg, paramt.mu(k,:), paramt.C(:,:,k));
    xxd = reshape(xxd, size(xx1));
    hold on
    contour(xx1, xx2, xxd, [1 1], '--k', 'linewidth', 2)
    hold off
end

for k=1:size(param1.mu,1)
    xxd = util_distM(xg, param1.mu(k,:), param1.C(:,:,k));
    xxd = reshape(xxd, size(xx1));
    hold on
    contour(xx1, xx2, xxd, [1 1], '-k', 'linewidth', 2)
    hold off
end

title('standard EM', 'FontSize', 15)
axis([min(paramt.mu(:,1))-1, tc+1, min(paramt.mu(:,2))-1, tc+1]), axis square


% truncated and censored data EM
subplot(122)
% data points
plot(x(:,1), x(:,2), 'xr', ...
    'markersize', 5)
hold on
plot(x(idxcx,1), x(idxcx,2), 'xr', ...
    'markersize', 5)
hold off

% centroids
hold on
plot(paramt.mu(:,1), paramt.mu(:,2), '+k', ...
    param2.mu(:,1), param2.mu(:,2), 'ok', 'linewidth', 2)
hold off

% contours
for k=1:size(paramt.mu,1)
    xxd = util_distM(xg, paramt.mu(k,:), paramt.C(:,:,k));
    xxd = reshape(xxd, size(xx1));
    hold on
    contour(xx1, xx2, xxd, [1 1], '--k', 'linewidth', 2)
    hold off
end

for k=1:size(param2.mu,1)
    xxd = util_distM(xg, param2.mu(k,:), param2.C(:,:,k));
    xxd = reshape(xxd, size(xx1));
    hold on
    contour(xx1, xx2, xxd, [1 1], '-k', 'linewidth', 2)
    hold off
end

title('truncated and censored EM', 'FontSize', 15)
axis([min(paramt.mu(:,1))-1, tc+1, min(paramt.mu(:,2))-1, tc+1]), axis square

