Num = 20; % the number of samples
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

X = xlsread('demo_data/X_sample.csv');
Yy_random = xlsread('demo_data/data_train.xlsx','Yy');
Mm_random = xlsread('demo_data/data_train.xlsx','M');
Yy1 = xlsread('demo_data/data_test.xlsx','Yy');
Mm1 = xlsread('demo_data/data_test.xlsx','M');

% predicting functional responses and calculating RMSE
Xx = X(1:Num,:);
Yy_random = Yy_random(1:60,1:Num);
Mm_random = Mm_random(1:60,1:Num);
Seq = 1:1:Num;
SE_RY = []; YP = [];

tic;
for i = 1:1:5 % testing sample number
        Sq = Seq;
        Sq(i) = [];
        Sq_test = i;
        Sq_train = Sq;
        [Yp,RY] = rmse_random(Xx,Yy_random,Mm_random,Yy1,Mm1,D,Sq_train,Sq_test);
        YP = [YP;Yp];
        SE_RY = [SE_RY;RY];
end
t = toc;
fprintf('the total running time is %7.4f s\n',t);

RMSE = sqrt(mean(SE_RY,2));

% plotting the predictive curve for testing sample

[Tp,Rp] = pol2cart(2*pi*Mm1(:,4),YP(4,:)');
[Tt,Rt] = pol2cart(2*pi*Mm1(:,4),Yy1(:,4));

plot([Tp;Tp(1)],[Rp;Rp(1)],'r');
hold on
plot(Tt,Rt,'*b');
