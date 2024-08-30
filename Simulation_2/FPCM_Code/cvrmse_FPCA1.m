% FPCM for closed curves in regular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train.csv');
Num = 20;
X = X(1:Num,:);
RMSE_FPCA = [];
for k = 1:100
    Yy_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');

    Yy_fixed = Yy_fixed(:,1:Num);
    Mm_fixed = Mm_fixed(:,1:Num);
    tic;
    SE_RY = rmse_FPCA1(X,X,Yy_fixed,Yy1,Mm_fixed,Mm1);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    % calculate cross-validated RMSE
    cv_mse = mean(SE_RY);
    cv_rmse = sqrt(cv_mse);
    RMSE_FPCA = [RMSE_FPCA;cv_rmse];
end

% RMSE = sqrt(mean(RMSE_FPCA.^2,2));
% boxplot(RMSE);
% ylim([0.25,0.75]);
% xlswrite('RMSE_FPCA_R_regular.xlsx',RMSE_FPCA);