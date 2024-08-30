
X = xlsread('demo_data\X_sample.xlsx');
Yy_random = xlsread('demo_data\data_random_train.xlsx','Yy');
Mm_random = xlsread('demo_data\data_random_train.xlsx','M');
Yy1 = xlsread('demo_data\data_test.xlsx','Yy');
Mm1 = xlsread('demo_data\data_test.xlsx','M');

% predicting functional responses and calculating RMSE
Num = 20; % the number of samples
Xx = X(1:Num,:);
Yy_random = Yy_random(:,1:Num); 
Mm_random = Mm_random(:,1:Num);
Seq = 1:1:Num;
SE_RY = []; YP = [];

tic;
for i = 1:1:5
        Sq = Seq;
        Sq(i) = [];
        Sq_test = i;
        Sq_train = Sq;
        [Yp,RY] = rmse_random(Xx,Yy_random,Mm_random,Yy1,Mm1,Sq_train,Sq_test);
        YP = [YP;Yp];
        SE_RY = [SE_RY;RY];
end
t = toc;
fprintf('the total running time for the 5 iterations is %7.4f s\n',t);

RMSE_RAND = sqrt(mean(SE_RY,2));

% plotting the predictive curve for testing sample
plot(Mm1(:,2)',YP(2,:),'r');
hold on
plot(Mm1(:,2),Yy1(:,2),'*b');
