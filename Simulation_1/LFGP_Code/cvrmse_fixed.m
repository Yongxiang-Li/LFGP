addpath('../simulation data/');

X = xlsread('../simulation data/X1_train_general1.xlsx');
Num = 20;


RMSE_FIX = [];
for k = 1:100
    Yy_fixed = xlsread(['../simulation data/data_fixed_rmse_general1\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_fixed = xlsread(['../simulation data/data_fixed_rmse_general1\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Xx = X(1:Num,:);
    Yy_fixed = Yy_fixed(:,1:Num); Mm_fixed = Mm_fixed(:,1:Num);
    Yy1 = Yy1(:,1:Num); Mm1 = Mm1(:,1:Num);

    % calculate cross-validated RMSE
    Seq = 1:1:Num; % sample sequence
    SE_FY = [];
    tic;
    for i = 1:1:Num
            Sq = Seq;
            Sq(i) = [];
            Sq_test = i;
            Sq_train = Sq;
            FY = rmse_fixed(Xx,Yy_fixed,Mm_fixed,Yy1,Mm1,Sq_train,Sq_test);
            SE_FY = [SE_FY;FY];
    end
    t = toc;
    fprintf('the total running time for the 20 iterations is %7.4f s\n',t);
    cv_mse1 = mean(SE_FY);
    cv_rmse1 = sqrt(cv_mse1);
    RMSE_FIX = [RMSE_FIX;cv_rmse1];
end
% xlswrite('RMSE_FIX_40_R_2_by.xlsx',RMSE_FIX);
% boxplot(RMSE_FIX);