addpath('../Data_total/');

X = xlsread('../Data_total/Input_Actuator.xlsx');
Xt = xlsread('../Data_total/Input_Actuator_test.xlsx');

Yy = xlsread('../Data_total/deformation_total.xlsx','Y');
M = xlsread('../Data_total/deformation_total.xlsx','M');

Yt = xlsread('../Data_total/deformation_total_t.xlsx','Y');
Mt = xlsread('../Data_total/deformation_total_t.xlsx','M');

Yy = 100*Yy; % unit conversion
Yt = 100*Yt; 

[SE_RY,YP] = rmse_FPCA_1(X,Xt,Yy,Yt,M,Mt);
% calculate cross-validated RMSE
cv_mse = mean(SE_RY,2);
cv_rmse = sqrt(cv_mse);
RMSE_FPCA = cv_rmse';

% boxplot(RMSE_FPCA);
% ylim([0.2,2]);
% xlswrite('RMSE_FPCA_TOTAL.xlsx',RMSE_FPCA);
% xlswrite('TP_FPCA.xlsx',YP);