
%clearvars -except aa smooth_cri
load('rawUv_5steps.mat');
[a,b_train]=max(v_train');   %v_train is the 32by1 vector
b_train=b_train';   %here b is the training output
[a,b_validation]=max(v_validation');
b_validation=b_validation';

mdlTreeBag = TreeBagger(100, U_train, b_train, 'method', 'classification', ...
                       'oobpred', 'on', 'minleaf', 30);

predict_v = predict(mdlTreeBag, U_test);
for k=1:length(predict_v)
    array_val(k,1)=str2num(predict_v{k});
end
array_val=array_val*0.25-0.125;
err=v_test-array_val;
rmse=sqrt(sum(err.^2)/length(v_test));

smooth_cri=sum(abs(diff(array_val)))/length(diff(array_val));

