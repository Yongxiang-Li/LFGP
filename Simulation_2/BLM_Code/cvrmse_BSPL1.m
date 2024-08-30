% BLM for closed curves in regular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train.csv');
Num = 20;
X = X(1:Num,:);
RMSE_BSPL = [];
for k = 1:100
    Yy_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'M');
    Yy = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Yy_fixed = Yy_fixed(:,1:Num); Mm_fixed = Mm_fixed(:,1:Num);
    Yy = Yy(:,1:Num); Mm = Mm(:,1:Num);
    tic;
    SE_RY = rmse_BSPL1(X,Yy_fixed,Mm_fixed,Yy,Mm);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    % calculate cross-validated RMSE
    cv_mse = mean(SE_RY);
    cv_rmse = sqrt(cv_mse);
    RMSE_BSPL = [RMSE_BSPL;cv_rmse];
end

% RMSE = sqrt(mean(RMSE_BSPL.^2,2));
% boxplot(RMSE);
% ylim([0.25,0.7]);
% xlswrite('RMSE_BSPL_R_regular.xlsx',RMSE_BSPL);