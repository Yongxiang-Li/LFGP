% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for closed curves in irregular grid condition

RMSE_RAND = xlsread('LFGP_CLOSED_Code\RMSE_RAND_CLOSED.xlsx');
RMSE_RAND_id = xlsread('LFGP_CLOSED_id_Code\RMSE_RAND_CLOSED.xlsx');
RMSE_FPCA = xlsread('FPCM_Code\RMSE_FPCA_R.xlsx');
RMSE_KRIG = xlsread('BK_Code\RMSE_KRIG_R.xlsx');
RMSE_BSPL = xlsread('BLM_Code\RMSE_BSPL_R.xlsx');

Seq = 1:51;

MSE_RAND = mean(RMSE_RAND(:,Seq).^2,2);
MSE_RAND_id = mean(RMSE_RAND_id(:,Seq).^2,2);
MSE_FPCA = mean(RMSE_FPCA.^2,2);
MSE_KRIG = mean(RMSE_KRIG.^2,2);
MSE_BSPL = mean(RMSE_BSPL.^2,2);

MSE = [MSE_RAND,MSE_RAND_id,MSE_BSPL,MSE_FPCA,MSE_KRIG];
RMSE = sqrt(MSE);

figure
boxplot(RMSE, 'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylim([0.25,0.7]);
ylabel('RMSE','fontsize',16);

% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for closed curves in regular grid condition

RMSE_FIX = xlsread('LFGP_CLOSED_Code\RMSE_FIX_CLOSED.xlsx');
RMSE_FIX_id = xlsread('LFGP_CLOSED_id_Code\RMSE_FIX_CLOSED.xlsx');
RMSE_FPCA = xlsread('FPCM_Code\RMSE_FPCA_R_regular.xlsx');
RMSE_KRIG = xlsread('BK_Code\RMSE_KRIG_R_regular.xlsx');
RMSE_BSPL = xlsread('BLM_Code\RMSE_BSPL_R_regular.xlsx');

Seq = 1:51;

MSE_FIX = mean(RMSE_FIX(:,Seq).^2,2);
MSE_FIX_id = mean(RMSE_FIX_id(:,Seq).^2,2);
MSE_FPCA = mean(RMSE_FPCA.^2,2);
MSE_KRIG = mean(RMSE_KRIG.^2,2);
MSE_BSPL = mean(RMSE_BSPL.^2,2);

MSE1 = [MSE_FIX,MSE_FIX_id,MSE_BSPL,MSE_FPCA,MSE_KRIG];
RMSE1 = sqrt(MSE1);

figure
boxplot(RMSE1, 'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylim([0.25,0.75]);
ylabel('RMSE','fontsize',16);