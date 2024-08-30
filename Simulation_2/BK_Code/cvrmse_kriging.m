% BK for closed curves in irregular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train.csv');
RMSE = [];
for k = 1:100
    M = xlsread(['../simulation data\data_random_rmse\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'M');
    Y = xlsread(['../simulation data\data_random_rmse\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'Yy');
    Mt = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Yt = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    tic;
    rmse = blind_kriging(X, M, Y, Mt, Yt);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    RMSE = [RMSE;rmse];
end

% xlswrite('RMSE_KRIG_R.xlsx',RMSE);