load PedDataScreen
global aic bic
%PedDataScreen  = PedDataScreen(PedDataScreen.Speed>1,:);
testlb = 1;
testub = 10;
PedDataScreen1=PedDataScreen(randperm(size(PedDataScreen,1)),:);%random rows for train
TRAIN=PedDataScreen(10001:43381,:);
TEST=PedDataScreen(1:10000,:);
PedTrainX=TRAIN(:,[8,10,19,27]);
PedTrainY=TRAIN(:,17);
PedTestX=TEST(:,[8,10,19,27]);
PedTestY=TEST(:,17);

% 
PedTestX = table2array(PedTestX);
PedTestY = table2array(PedTestY);
PedTrainX = table2array(PedTrainX);
PedTrainY = table2array(PedTrainY);
y = [PedDataScreen.Range,PedDataScreen.Transversal,PedDataScreen.Speed,PedDataScreen.vp,PedDataScreen.Ax];
K=13;
options = statset('MaxIter',2000);
res_obj = fitgmdist(y, K, 'Options',options);
aic(K) = res_obj.AIC;
bic(K) = res_obj.BIC;
param.pp = res_obj.PComponents;
param.mu = res_obj.mu;
param.C = res_obj.Sigma;

pdf=0;

apregmm = zeros(testub-testlb+1,1);
apregame = zeros(testub-testlb+1,1);
for i = testlb:testub
    x1 = PedTestX(i,1);
    x2 = PedTestX(i,2);
    x3 = PedTestX(i,3);
    x4 = PedTestX(i,4);
%    [x1,x2,x3,x4] = PedTestX(i,:);
    cond = 0;
    for j = 1:K
        gain = param.pp(j)*mvnpdf([x1,x2,x3,x4],param.mu(j,1:4),param.C(1:4,1:4,j));
        cond = cond + gain;
    end
    a5 = zeros(501,1);
    for ii = 1:501
        x5 = ii/100-3.01;
        for j = 1:K
            Vq = param.pp(j)*mvnpdf([x1,x2,x3,x4,x5],param.mu(j,:),param.C(:,:,j));
            a5(ii) = a5(ii)+Vq;
        end
        a5(ii) = a5(ii)/cond;
    end
    [~,apregmm(i-testlb+1)] = max(a5);
    apregmm(i-testlb+1) = apregmm(i-testlb+1)/100-3.01;
%    apregame(i-testlb+1) = 2*x1/16.5776 - x3;
    %apregame(i-testlb+1) = -0.0151*x1+0.0390*x2-0.027*x3+0.03*x4+0.2469*x3^2/x1-0.0116*x4*x2;
    apregame(i-testlb+1) = 2*0.0272*-0.0594*x1*x4/(x4+2*0.0594*x2)-0.0272*x3;
%     [-0.0594555520593611,0.0271862198328049]
%     a = 2*alpha(2)*alpha(1).*ux(:,1).*ux(:,4)./(ux(:,4)-2*alpha(1).*ux(:,2))-alpha(2)*ux(:,3);
end


%mdlTree = RegressionTree.fit(PedTrainX,PedTrainY,'PredictorNames',labels);
mdlTree = RegressionTree.fit(PedTrainX,PedTrainY);
forecast(1).Y = predict(mdlTree,PedTestX);
net = fitnet(20);
net = train(net, PedTrainX', PedTrainY');
forecast(2).Y = net(PedTestX')';
mdlTreeBag = TreeBagger(100, PedTrainX, PedTrainY, 'method', 'regression', ...
                       'oobpred', 'on', 'minleaf', 30);
forecast(3).Y = predict(mdlTreeBag, PedTestX);

%idx = testDates > datenum('Jun-01-2008') & testDates < datenum('Jul-01-2008');

%Dates = testDates(idx);
idx = testlb:testub;

figure('Units','Normalized','Position',[0.05,0.4,0.4,0.5]) %subplot(2,1,1)
% 
% hPlot1 = plot(idx, [PedTestY(idx),forecast(1).Y(idx),...
%     forecast(2).Y(idx),forecast(3).Y(idx)],'LineWidth',2);
% set(hPlot1(1),'LineWidth',5,'Color',[1 1 0],'DisplayName','Actual');
% set(hPlot1(2),'Color',[0 1 0],'DisplayName','Regression Tree');
% set(hPlot1(3),'DisplayName','Neural Network');
% set(hPlot1(4),'Color',[0.5 0.5 0.5],'DisplayName','Bagged Regression Trees');
% legend('show'),
% %datetick('x','mmm-dd','keepticks'), xlabel('Time'),
% ylabel('Load'),
% title('Ax Prediction','FontSize',12,'FontWeight','Bold')
% 
% subplot(2,1,2)
% hPlot2 = plot(idx,[PedTestY(idx)-forecast(1).Y(idx)...
%     PedTestY(idx)-forecast(2).Y(idx),PedTestY(idx)-forecast(3).Y(idx)]);
% set(hPlot2(1),'Color',[0 1 0],'DisplayName','Regression Tree Error');
% set(hPlot2(2),'Color','r','DisplayName','Neural Network Error');
% set(hPlot2(3),'Color',[0.5 0.5 0.5],'DisplayName','Bagged Regression Trees Error');
% %datetick('x','mmm-dd','keepticks'), xlabel('Time'),
% title('Prediction Error','FontSize',12,'FontWeight','Bold')
% ylabel('Residuals'), grid on
% legend('show')


hPlot1 = plot(idx, [PedTestY(idx),forecast(2).Y(idx),...
    forecast(3).Y(idx),apregmm,apregame],'LineWidth',2);
set(hPlot1(1),'LineWidth',5,'Color',[1 1 0],'DisplayName','Actual');

set(hPlot1(2),'DisplayName','Neural Network');
set(hPlot1(3),'Color',[0.5 0.5 0.5],'DisplayName','Bagged Regression Trees');
set(hPlot1(4),'Color',[0 0 1],'DisplayName','GMM');
set(hPlot1(5),'Color',[1 0 1],'DisplayName','Nash Equilibrium');
legend('show'),
%datetick('x','mmm-dd','keepticks'), xlabel('Time'),
ylabel('Ax'),
title('Ax Prediction','FontSize',12,'FontWeight','Bold')

% subplot(2,1,2)
% hPlot2 = plot(idx,[...
%     PedTestY(idx)-forecast(2).Y(idx),PedTestY(idx)-forecast(3).Y(idx)]);
% 
% set(hPlot2(1),'Color','r','DisplayName','Neural Network Error');
% set(hPlot2(2),'Color',[0.5 0.5 0.5],'DisplayName','Bagged Regression Trees Error');
% %datetick('x','mmm-dd','keepticks'), xlabel('Time'),
% title('Prediction Error','FontSize',12,'FontWeight','Bold')
% ylabel('Residuals'), grid on
% legend('show')
disp('mean:')
disp(mean(PedTestY(idx)-forecast(2).Y(idx)));
disp(mean(PedTestY(idx)-forecast(3).Y(idx)));
disp('R')
disp(corrcoef(PedTestY(idx),forecast(2).Y(idx)));
disp(corrcoef(PedTestY(idx),forecast(3).Y(idx)));
disp('std:');
disp(std(PedTestY(idx)-forecast(2).Y(idx)));
disp(std(PedTestY(idx)-forecast(3).Y(idx)));

sqrt(sum((PedTestY(idx)-forecast(2).Y(idx)).^2)/length(idx))
sqrt(sum((PedTestY(idx)-forecast(3).Y(idx)).^2)/length(idx))
sqrt(sum((PedTestY(idx)-apregmm).^2)/length(idx))
sqrt(sum((PedTestY(idx)-apregame).^2)/length(idx))