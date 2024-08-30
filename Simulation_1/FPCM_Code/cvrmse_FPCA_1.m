% FPCM for functional outputs under irregular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train_general1.xlsx');
Num = 20;
X = X(1:Num,:);
RMSE_FPCA = [];
for k = 1:100
    Yy_random = xlsread(['../simulation data/data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_random = xlsread(['../simulation data/data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');

    Yy_random = Yy_random(:,1:Num);
    Mm_random = Mm_random(:,1:Num);
    tic;
    SE_RY = rmse_FPCA_1(X,X,Yy_random,Yy1,Mm_random,Mm1);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    % calculate cross-validated RMSE
    cv_mse = mean(SE_RY);
    cv_rmse = sqrt(cv_mse);
    RMSE_FPCA = [RMSE_FPCA;cv_rmse];
end

% RMSE = sqrt(mean(RMSE_FPCA.^2,2));
% boxplot(RMSE); 
% ylim([0.37,0.48]);
% xlswrite('RMSE_FPCA_R_1.xlsx',RMSE_FPCA);