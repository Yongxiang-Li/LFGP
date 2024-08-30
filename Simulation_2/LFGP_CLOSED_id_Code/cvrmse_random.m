addpath('../simulation data/');

X = xlsread('../simulation data\X1_train.csv');
Num = 20;
p = 20; % the number of control points
d = 3; % the order of B-spline
knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';
transfer_matrix = coeff_transform_matrix(knots, p, d);
bse = basismatrix(d, p, knots, [0,1]);
b = bse(1,:)-bse(2,:);
dd = d-1;    dp = p-1;    dknots = knots(2:end-1);
dbse = basismatrix(dd, dp, dknots, [0,1]);
db = dbse(1,:)-dbse(2,:);
db = db*transfer_matrix;
st = [b;db];
D = st; % constraint matrix D

RMSE_RAND = [];
for k = 1:100
    Yy_random = xlsread(['../simulation data\data_random_rmse\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_random = xlsread(['../simulation data\data_random_rmse\simulation_data_constraint_random_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');

    % calculate cross-validated RMSE
    Xx = X(1:Num,:);
    Yy_random = Yy_random(1:60,1:Num);
    Mm_random = Mm_random(1:60,1:Num);
    Seq = 1:1:Num; % sample sequence
    SE_RY = [];
    tic;
    for i = 1:1:Num
            Sq = Seq;
            Sq(i) = [];
            Sq_test = i;
            Sq_train = Sq;
            RY = rmse_random(Xx,Yy_random,Mm_random,Yy1,Mm1,D,Sq_train,Sq_test);
            SE_RY = [SE_RY;RY];
    end
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    SE = sqrt(mean(SE_RY,2));
    ind = find(SE<=0.6);
    cv_mse = mean(SE_RY(ind',:),1);
    cv_rmse = sqrt(cv_mse);
    RMSE_RAND = [RMSE_RAND;cv_rmse];
end
% xlswrite('RMSE_RAND_CLOSED.xlsx',RMSE_RAND);
% boxplot(RMSE_RAND);