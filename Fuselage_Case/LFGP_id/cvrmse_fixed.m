addpath('../Data_total/');

X = xlsread('../Data_total/Input_Actuator.xlsx');
Xt = xlsread('../Data_total/Input_Actuator_test.xlsx');
X = -X; Xt = -Xt;
Num = 25; % the number of test samples
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


Yy_random = xlsread('../Data_total/deformation_total_fixed.xlsx','Y');
Mm_random = xlsread('../Data_total/deformation_total_fixed.xlsx','M');

Yt_random = xlsread('../Data_total/deformation_total_t.xlsx','Y');
Mt_random = xlsread('../Data_total/deformation_total_t.xlsx','M');

Yy_random = 100*Yy_random; % unit conversion
Yt_random = 100*Yt_random; 


% calculate cross-validated RMSE

Seq = 1:1:30; % training sample sequence
SE_RY = [];
YP = [];
PARA = [];

for i = 1:1:Num
        Sq_test = i;
        Sq_train = Seq;
        [RY,Yp,model] = rmse_fixed(X,Xt,Yy_random,Yt_random,Mm_random,Mt_random,D,Sq_train,Sq_test);
        YP = [YP,Yp'];
        SE_RY = [SE_RY;RY];
        para = [model.theta,model.delta];
        PARA = [PARA;para];
end

cv_mse = mean(SE_RY,2);
cv_rmse = sqrt(cv_mse);
RMSE_FIX = cv_rmse';

% boxplot(RMSE_FIX);
% ylim([0.2,2.4]);
% xlswrite('RMSE_LFGP_TOTAL_regular.xlsx',RMSE_FIX);
% xlswrite('TP_LFGP.xlsx',YP);