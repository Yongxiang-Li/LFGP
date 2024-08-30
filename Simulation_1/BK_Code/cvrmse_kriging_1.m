% BK for functional outputs under irregular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train_general1.xlsx');
RMSE = [];
for k = 1:100
    M = xlsread(['../simulation data\data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'M');
    Y = xlsread(['../simulation data\data_random_rmse_general1\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'Yy');
    Mt = xlsread(['../simulation data\data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Yt = xlsread(['../simulation data\data_test_rmse_general1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    tic;
    rmse = blind_kriging1(X, M, Y, Mt, Yt);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    RMSE = [RMSE;rmse];
end

% xlswrite('RMSE_KRIG_R_1_by.xlsx',RMSE);