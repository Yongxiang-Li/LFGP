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

RMSE_FIX = [];
for k = 1:100
    Yy_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'Yy');
    Mm_fixed = xlsread(['../simulation data\data_fixed_rmse\simulation_data_constraint_fixed_rmse',num2str(k),'.xlsx'],'M');
    Yy1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'Yy');
    Mm1 = xlsread(['../simulation data\data_test_rmse1\simulation_data_constraint_test_rmse',num2str(k),'.xlsx'],'M');

    % calculate cross-validated RMSE
    Xx = X(1:Num,:);
    Yy_fixed = Yy_fixed(:,1:Num);
    Mm_fixed = Mm_fixed(:,1:Num);
    Seq = 1:1:Num; % sample sequence
    SE_FY = [];
    tic;
    for i = 1:1:Num
            Sq = Seq;
            Sq(i) = [];
            Sq_test = i;
            Sq_train = Sq;
            RY = rmse_fixed(Xx,Yy_fixed,Mm_fixed,Yy1,Mm1,D,Sq_train,Sq_test);
            SE_FY = [SE_FY;RY];
    end
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    SE = sqrt(mean(SE_FY,2));
    ind = find(SE<=0.75);
    cv_mse = mean(SE_FY(ind',:),1);
    cv_rmse = sqrt(cv_mse);
    RMSE_FIX = [RMSE_FIX;cv_rmse];
end
% xlswrite('RMSE_FIX_CLOSED.xlsx',RMSE_FIX);
% boxplot(RMSE_FIX);