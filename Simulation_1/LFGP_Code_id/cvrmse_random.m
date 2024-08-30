addpath('../simulation data/');

X = xlsread('../simulation data/X1_train_general1.xlsx');
Num = 20;


RMSE_RAND = [];
for k = 1:100
    Yy_random = xlsread(['../simulation data/data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_random = xlsread(['../simulation data/data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data/data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Xx = X(1:Num,:);
    Yy_random = Yy_random(:,1:Num); Mm_random = Mm_random(:,1:Num);
    Yy1 = Yy1(:,1:Num); Mm1 = Mm1(:,1:Num);

    % calculate cross-validated RMSE
    Seq = 1:1:Num; % sample sequence
    SE_RY = [];
    tic;
    for i = 1:1:Num
            Sq = Seq;
            Sq(i) = [];
            Sq_test = i;
            Sq_train = Sq;
            RY = rmse_random(Xx,Yy_random,Mm_random,Yy1,Mm1,Sq_train,Sq_test);
            SE_RY = [SE_RY;RY];
    end
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    cv_mse = mean(SE_RY);
    cv_rmse = sqrt(cv_mse);
    RMSE_RAND = [RMSE_RAND;cv_rmse];
end
% xlswrite('RMSE_RAND_40_R_2_by.xlsx',RMSE_RAND);
% boxplot(RMSE_RAND);
% RMSE = sqrt(mean(RMSE_RAND.^2,2));
% boxplot(RMSE);
% ylim([0.37,0.48]);