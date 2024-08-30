addpath('../Data_total/');

X = xlsread('../Data_total/Input_Actuator.xlsx');
Xt = xlsread('../Data_total/Input_Actuator_test.xlsx');
X = -X; Xt = -Xt;

Yy = xlsread('../Data_total/deformation_total_fixed.xlsx','Y');
Mm = xlsread('../Data_total/deformation_total_fixed.xlsx','M');

Yt = xlsread('../Data_total/deformation_total_t.xlsx','Y');
Mt = xlsread('../Data_total/deformation_total_t.xlsx','M');

Yy = 100*Yy; % unit conversion
Yt = 100*Yt; 

[SE_RY,YP] = rmse_BSPL_2(X,Xt,Yy,Yt,Mm,Mt);
% calculate cross-validated RMSE
cv_mse = mean(SE_RY,2);
cv_rmse = sqrt(cv_mse);
RMSE_BSPL = cv_rmse';

% boxplot(RMSE_BSPL);
% ylim([0.2,2.4]);
% xlswrite('RMSE_BSPL_TOTAL_regular.xlsx',RMSE_BSPL);
% xlswrite('TP_BSPL.xlsx',YP);