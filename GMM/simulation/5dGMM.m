% data
clear
load PedDataScreen
global aic bic
%PedDataScreen = PedDataScreen(PedDataScreen.Speed > 0.1,:);

y = [PedDataScreen.Range,PedDataScreen.Transversal,PedDataScreen.Speed,PedDataScreen.vp];

%'XBinLimits',[0 1]

% truncation and censoring
% tc = Inf;
% idxt = y(:,1)>0;    % truncation
% idxcu = y>tc;       % censoring
% idxcl = y<-Inf;        % censoring
% idxc = any(idxcu | idxcl,2); %censoring two way
% 
% x = y;
% x(idxcu) = tc;
% x(idxcl) = -Inf;
% x = x(idxt,:);
% idxcx = idxc(idxt);
% 
% pattern = em_censor_pattern(x, [-Inf; tc]);
% xc = x(pattern.complete,:);

% 
% model fitting standard
% for K = 1:12
% 
% % initialize by k-means 
% 
K=13;
% par_k = init_kmeans(x, K, 5);
% par_k_std_ll = zeros(5,1);
% for m=1:5
%     [temp,par_k_std_ll(m)] = em_std_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C);
% end
% [temp,par_k_std_idx] = max(par_k_std_ll);
% param1_init = par_k{par_k_std_idx};
% 
% param1 = em_std(x, K, param1_init);
options = statset('MaxIter',2000);
%obj = gmdistribution.fit(X,2,'Options',options);
%x=y(:,1:2);
%res_obj = fitgmdist(x, K, 'Start', str_start,'Options',options);
res_obj = fitgmdist(y, K, 'Options',options);
aic(K) = res_obj.AIC;
bic(K) = res_obj.BIC;
param.pp = res_obj.PComponents;
param.mu = res_obj.mu;
param.C = res_obj.Sigma;




% figure;  
% y = [zeros(10000,1);ones(10000,1)];  
% hold on;  
% ezsurf (@(x1,x2)pdf(res_obj,[x1 x2]),get(gca,{'XLim','YLim'}));  
% title('{\bf 3-D Fitted Gaussian Mixture}');  
% hold off;  

% end
% 
% K = 8;
% par_k = init_kmeans(x, K, 5);
% 
% par_k_tc_ll = zeros(5,1);
% for m=1:5
%     [temp,par_k_tc_ll(m)] = em_tc1_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C, pattern, [0;Inf]);
% end
% [temp,par_k_tc_idx] = max(par_k_tc_ll);
% param2_init = par_k{par_k_tc_idx};
% param2 = em_tc1(x, K, [0;Inf], [-Inf; tc], param2_init);
% [temp,idx] = sort(param2.mu(:,1), 'ascend');
% param2.pp = param2.pp(idx);
% param2.mu = param2.mu(idx,:);
% param2.C = param2.C(:,:,idx);






% plot(bic);
% for i = 2:length(bic)
%     diff(i) = (bic(i) - bic(i-1))/abs(bic(i-1));
% end
% plot(diff)
% standard EM

%[X,Y]=meshgrid(linspace(-2,10,100),linspace(0,0.3,100));


% LEFT = 1;
% for i = 1:K
%     trunc(i) = cdf('Normal',0,param1.mu(i,1),sqrt(param1.C(1,1,i)));
%     LEFT = LEFT - param2.pp(i)*trunc(i);
% end

% for i = 1:K
% %     d(i) = param2.pp(i)*(1-trunc(i))/LEFT;
% %     multi(i) = 1/(1-trunc(i))*d(i);
%     p{i}=mvnpdf([X(:) Y(:)],param.mu(i,1:2),param.C(1:2,1:2,i));%????????????Z?
%     p{i}=reshape(p{i},size(X));%?Z??????????
%     p{i} = p{i}.*param.pp(i);
%     pdf=pdf+p{i};
% end
% surf(X,Y,pdf,'EdgeColor','none')
%RANGE TRANSVERSAL SPEED VP AX
pdf=0;
%[X1,X2,X3,X4,X5]=ndgrid(linspace(0,30,20),linspace(-10,10,20),linspace(0,10,20),linspace(-4,4,20),linspace(-2,2,20));


% for i = 1:K
%     trunc(i) = 1-(1-cdf('Normal',0,param1.mu(i,1),sqrt(param1.C(1,1,i)))) * ...
%         (1-cdf('Normal',0,param1.mu(i,2),sqrt(param1.C(2,2,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,3),sqrt(param1.C(3,3,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,4),sqrt(param1.C(4,4,i))));
%     LEFT = LEFT - param1.pp(i)*trunc(i);
% end

% for i = 1:K
%     d(i) = param.pp(i);
%     p{i}=mvnpdf([X1(:),X2(:),X3(:),X4(:),X5(:)],param.mu(i,:),param.C(:,:,i));
%     p{i}=reshape(p{i},size(X1));
%     p{i} = p{i}.*d(i);
%     pdf=pdf+p{i};
% end
x1 = 1;
x2 = 0;
x3 = 3;

cond = 0;
for j = 1:K
    gain = param.pp(j)*mvnpdf([x1,x2,x3],param.mu(j,1:3),param.C(1:3,1:3,j));
    cond = cond + gain;
end

a5 = zeros(501,1);
%[x1,x2,x3,x4] = [10,2,3,2];
for i = 1:501
    x4 = i/100-3.01;
%    [x1,x2,x3,x4] = ndgrid(0:0.01:2, 0.07, 1.20, 0.05);
%     [x2,x3,x4] = ndgrid(0.13, 1.20, 0.05); 
%     Vq = interpn(X1,X2,X3,X4,pdf,x1,x2,x3,x4);
    for j = 1:K
        Vq = param.pp(j)*mvnpdf([x1,x2,x3,x4],param.mu(j,:),param.C(:,:,j));
        a5(i) = a5(i)+Vq;
    end
    a5(i) = a5(i)/cond;
end
plot(-3:0.01:2,a5)
% hold on
% 
% x1 = 10;
% x2 = -2;
% x3 = 2;
% x4 = 1;
% cond = 0;
% for j = 1:K
%     gain = param.pp(j)*mvnpdf([x1,x2,x3,x4],param.mu(j,1:4),param.C(1:4,1:4,j));
%     cond = cond + gain;
% end
% 
% a5 = zeros(501,1);
% for i = 1:501
%     x5 = i/100-3.01;
%     for j = 1:K
%         Vq = param.pp(j)*mvnpdf([x1,x2,x3,x4,x5],param.mu(j,:),param.C(:,:,j));
%         a5(i) = a5(i)+Vq;
%     end
%     a5(i) = a5(i)/cond;
% end
% plot(-3:0.01:2,a5)
% hold on
% 
% x1 = 15;
% x2 = -2;
% x3 = 4;
% x4 = 1;
% cond = 0;
% for j = 1:K
%     gain = param.pp(j)*mvnpdf([x1,x2,x3,x4],param.mu(j,1:4),param.C(1:4,1:4,j));
%     cond = cond + gain;
% end
% 
% a5 = zeros(501,1);
% for i = 1:501
%     x5 = i/100-3.01;
%     for j = 1:K
%         Vq = param.pp(j)*mvnpdf([x1,x2,x3,x4,x5],param.mu(j,:),param.C(:,:,j));
%         a5(i) = a5(i)+Vq;
%     end
%     a5(i) = a5(i)/cond;
% end
% plot(-3:0.01:2,a5)
% hold on



% 
% pdf3=0;
% LEFT = 1;
% for i = 1:K
%     trunc(i) = 1-(1-cdf('Normal',0,param1.mu(i,1),sqrt(param1.C(1,1,i)))) * ...
%         (1-cdf('Normal',0,param1.mu(i,2),sqrt(param1.C(2,2,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,3),sqrt(param1.C(3,3,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,4),sqrt(param1.C(4,4,i))));
%     LEFT = LEFT - param1.pp(i)*trunc(i);
% end
% 
% for i = 1:K
%     d(i) = param1.pp(i)*(1-trunc(i))/LEFT;
%     multi(i) = 1/(1-trunc(i))*d(i);
%     p{i}=mvnpdf([X1(:),X2(:),X4(:)],param1.mu(i,[1:2,4]),param1.C([1:2,4],[1:2,4],i));
%     p{i}=reshape(p{i},size(X1));
%     p{i} = p{i}.*multi(i);
%     pdf3=pdf3+p{i};
% end


% [X2,X3,X4]=ndgrid(linspace(0,0.3,100),linspace(0,6,100),linspace(0,0.5,100));
% 
% pdf234=0;
% LEFT = 1;
% for i = 1:K
%     trunc(i) = 1-(1-cdf('Normal',0,param1.mu(i,1),sqrt(param1.C(1,1,i)))) * ...
%         (1-cdf('Normal',0,param1.mu(i,2),sqrt(param1.C(2,2,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,3),sqrt(param1.C(3,3,i))))* ...
%         (1-cdf('Normal',0,param1.mu(i,4),sqrt(param1.C(4,4,i))));
%     LEFT = LEFT - param1.pp(i)*trunc(i);
% end
% 
% for i = 1:K
%     d(i) = param1.pp(i)*(1-trunc(i))/LEFT;
%     multi(i) = 1/(1-trunc(i))*d(i);
%     p{i}=mvnpdf([X2(:),X3(:),X4(:)],param1.mu(i,2:4),param1.C(2:4,2:4,i));
%     p{i}=reshape(p{i},size(X2));
%     p{i} = p{i}.*multi(i);
%     pdf234=pdf234+p{i};
% end

% x2 : 0.07 0.1 0.13
% x3 : 1.2
% x4 : 0.1
% figure(3)
% for i = 1:201
%     x1 = i/100-0.01;
% %    [x1,x2,x3,x4] = ndgrid(0:0.01:2, 0.07, 1.20, 0.05);
%     [x2,x3,x4] = ndgrid(0.07, 1.20, 0.05); 
%     Vq = interpn(X1,X2,X3,X4,pdf,x1,x2,x3,x4);
%     v1(i) = Vq;
%     hold on
%     
% end
% plot(0:0.01:2,v1)
% hold on
% 
% for i = 1:201
%     x1 = i/100-0.01;
% %    [x1,x2,x3,x4] = ndgrid(0:0.01:2, 0.07, 1.20, 0.05);
%     [x2,x3,x4] = ndgrid(0.1, 1.20, 0.05); 
%     Vq = interpn(X1,X2,X3,X4,pdf,x1,x2,x3,x4);
%     v2(i) = Vq;
%     hold on
% end
% plot(0:0.01:2,v2)
% hold on
% 
% for i = 1:201
%     x1 = i/100-0.01;
% %    [x1,x2,x3,x4] = ndgrid(0:0.01:2, 0.07, 1.20, 0.05);
%     [x2,x3,x4] = ndgrid(0.13, 1.20, 0.05); 
%     Vq = interpn(X1,X2,X3,X4,pdf,x1,x2,x3,x4);
%     v3(i) = Vq;
%     hold on
% end
% plot(0:0.01:2,v3)
% hold on



