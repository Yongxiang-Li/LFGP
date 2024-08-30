% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for composite fuselages in irregular grid

RMSE_FOGP = xlsread('LFGP\RMSE_LFGP_TOTAL.xlsx');
RMSE_FOGP_id = xlsread('LFGP_id\RMSE_LFGP_TOTAL.xlsx');
RMSE_FPCA = xlsread('FPCM\RMSE_FPCA_TOTAL.xlsx');
RMSE_KRIG = xlsread('BK\RMSE_KRIG_TOTAL.xlsx');
RMSE_BSPL = xlsread('BLM\RMSE_BSPL_TOTAL.xlsx');

RMSE1 = [RMSE_FOGP',RMSE_FOGP_id',RMSE_BSPL',RMSE_FPCA',RMSE_KRIG'];
figure
boxplot(RMSE1,'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylabel('RMSE','fontsize',16);
ylim([0.2,2]);

% Predictive RMSE of LFGP, LFGP-IC, BLM, FPCM and BK for composite fuselages in regular grid

RMSE_FOGP = xlsread('LFGP\RMSE_LFGP_TOTAL_regular.xlsx');
RMSE_FOGP_id = xlsread('LFGP_id\RMSE_LFGP_TOTAL_regular.xlsx');
RMSE_FPCA = xlsread('FPCM\RMSE_FPCA_TOTAL_regular.xlsx');
RMSE_KRIG = xlsread('BK\RMSE_KRIG_TOTAL_regular.xlsx');
RMSE_BSPL = xlsread('BLM\RMSE_BSPL_TOTAL_regular.xlsx');

RMSE2 = [RMSE_FOGP',RMSE_FOGP_id',RMSE_BSPL',RMSE_FPCA',RMSE_KRIG'];
figure
boxplot(RMSE2,'Labels',{'LFGP','LFGP-IC','BLM','FPCM','BK'});
ylabel('RMSE','fontsize',16);
ylim([0.2,2.2]);

% Plot the predicted dimensions of composite fuselage

def_Yp = xlsread('LFGP\result\TP_LFGP.xlsx');
def_Yt = xlsread('Data_total\deformation_total_test.xlsx','Y');
Mt_def = xlsread('Data_total\deformation_total_test.xlsx','M');

Mt = Mt_def(:,7);  
Smp_def_Yt = 100*def_Yt(:,7)-190;
Smp_def_Yp = def_Yp(:,7)-190;

[Mt,J] = sort(Mt);  
Smp_def_Yt = Smp_def_Yt(J);
Smp_def_Yp = Smp_def_Yp(J);

Mm = 0:0.01:1; Mm = Mm(1,2:100);
Smp_def_Yt = interp1(Mt,Smp_def_Yt,Mm','spline');
Smp_def_Yp = interp1(Mt,Smp_def_Yp,Mm','spline');

location = xlsread('Data_total\Location.xlsx');
Y0 = location(:,3); Z0 = location(:,4);

[theta,rho] = cart2pol(Y0,Z0);
idx = find(theta>0);
idx1 = find(theta<0);
theta(idx) = 2*pi - theta(idx);
theta(idx1) = -theta(idx1);
[raio,I] = sort(theta);
Rho = rho(I); Rho = 100*Rho-190;

% transformation of coordinates
Ww = 2*pi*Mm';
[Y,Z] = pol2cart(raio,Rho);
[Yt,Zt] = pol2cart(Ww,Smp_def_Yt);
[Yp,Zp] = pol2cart(Ww,Smp_def_Yp);

Y = [Y;Y(1)]; Z = [Z;Z(1)]+1.5;
Yt = [Yt;Yt(1)]; Zt = [Zt;Zt(1)]+0.5;
Yp = [Yp;Yp(1)]; Zp = [Zp;Zp(1)];

figure
p1 = plot(Y,Z,'b--','LineWidth',1.4);
hold on
p2 = plot(Yt,Zt,'-.','Color',[0.2,0.8,0.2],'LineWidth',1.9);
hold on
p3 = plot(Yp,Zp,'r-','LineWidth',1.3);
hold off
axis equal
ylim([-70,90]);
xlabel('Y-direction','fontsize',12);
ylabel('Z-direction','fontsize',12);
legend([p1,p2,p3],{'Design Values','Actual Values','Predicted Values'},'location','southoutside','FontSize',10,'Orientation','horizontal');