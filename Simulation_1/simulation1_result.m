% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for general functional outputs in irregular grid condition

RMSE_RAND1 = xlsread('LFGP_Code\RMSE_RAND_40_R_2_by.xlsx');
RMSE_RAND_id = xlsread('LFGP_Code_id\RMSE_RAND_40_R_2_by.xlsx');
RMSE_FPCA1 = xlsread('FPCM_Code\RMSE_FPCA_R_1.xlsx');
RMSE_KRIG1 = xlsread('BK_Code\RMSE_KRIG_R_1_by.xlsx');
RMSE_BSPL1 = xlsread('BLM_Code\RMSE_BSPL_R_1.xlsx');

Seq = 1:51;

MSE_RAND1 = mean(RMSE_RAND1(:,Seq).^2,2);
MSE_RAND_id = mean(RMSE_RAND_id(:,Seq).^2,2);
MSE_FPCA1 = mean(RMSE_FPCA1.^2,2);
MSE_KRIG1 = mean(RMSE_KRIG1.^2,2);
MSE_BSPL1 = mean(RMSE_BSPL1.^2,2);

MSE1 = [MSE_RAND1,MSE_RAND_id,MSE_BSPL1,MSE_FPCA1,MSE_KRIG1];
RMSE1 = sqrt(MSE1);

figure
boxplot(RMSE1, 'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylim([0.37,0.47]);
ylabel('RMSE','fontsize',16);

% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for general functional outputs in regular grid condition

RMSE_FIX = xlsread('LFGP_Code\RMSE_FIX_40_R_2_by.xlsx');
RMSE_FIX_id = xlsread('LFGP_Code_id\RMSE_FIX_40_R_2_by.xlsx');
RMSE_FPCA_FIX = xlsread('FPCM_Code\RMSE_FPCA_R_1_regular.xlsx');
RMSE_KRIG_FIX = xlsread('BK_Code\RMSE_KRIG_R_1_regular.xlsx');
RMSE_BSPL_FIX = xlsread('BLM_Code\RMSE_BSPL_R_1_regular.xlsx');

Seq = 1:51;

MSE_FIX = mean(RMSE_FIX(:,Seq).^2,2);
MSE_FIX_id = mean(RMSE_FIX_id(:,Seq).^2,2);
MSE_FPCA_FIX = mean(RMSE_FPCA_FIX.^2,2);
MSE_KRIG_FIX = mean(RMSE_KRIG_FIX.^2,2);
MSE_BSPL_FIX = mean(RMSE_BSPL_FIX.^2,2);

MSE2 = [MSE_FIX,MSE_FIX_id,MSE_BSPL_FIX,MSE_FPCA_FIX,MSE_KRIG_FIX];
RMSE2 = sqrt(MSE2);

figure
boxplot(RMSE2, 'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylim([0.37,0.52]);
ylabel('RMSE','fontsize',16);