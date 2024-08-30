% BK for closed curves in regular grid conditions
addpath('../simulation data/');

X = xlsread('../simulation data\X1_train.csv');
RMSE = [];
for k = 1:100
    M = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'M');
    Y = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'Yy');
    Mt = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');
    Yt = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    M = M(21:40,:); Y = Y(21:40,:);
    tic;
    rmse = blind_kriging1(X, M, Y, Mt, Yt);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    RMSE = [RMSE;rmse];
end

% xlswrite('RMSE_KRIG_R_regular.xlsx',RMSE);