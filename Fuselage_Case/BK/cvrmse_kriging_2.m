% BK for fuselage shape under regular grid conditions
addpath('../Data_total/');

X = xlsread('../Data_total/Input_Actuator.xlsx');
Xt = xlsread('../Data_total/Input_Actuator_test.xlsx');
X = -X; Xt = -Xt;

M = xlsread('../Data_total/deformation_total_fixed.xlsx','M');
Y = xlsread('../Data_total/deformation_total_fixed.xlsx','Y');

Mt = xlsread('../Data_total/deformation_total_t.xlsx','M');
Yt = xlsread('../Data_total/deformation_total_t.xlsx','Y');

Y = 100*Y; % unit conversion
Yt = 100*Yt;

rmse = blind_kriging2(X, Xt, M, Y, Mt, Yt);
RMSE = rmse';

% boxplot(RMSE);
% ylim([0.2,2.4]);
% xlswrite('RMSE_KRIG_TOTAL_regular.xlsx',RMSE);
% xlswrite('TP_KRIG.xlsx',YP);