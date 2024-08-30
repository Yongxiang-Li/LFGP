function SE = rmse_FPCA(X,Xt,A,At,M,Mt)

A = A'; At = At';
b=size(A,2);
h = size(X,1);

% Normalization of data: obtain the normalized matrix SA
% MEAN = mean(A); STD = std(A);
for i=1:b
    SA(:,i)=(A(:,i)-mean(A(:,i)))/std(A(:,i));
end
MEANt = mean(At); STDt = std(At);

p = 10; % % the number of control points
d = 3; % the order of B-spline
knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';

C = [];
for j = 1:h
    ss = M(:,j);
    Y = SA(j,:);
    U = basismatrix(d, p, knots, ss');
    P = U\Y';
    C = [C,P];
end
C = C';
w = 0:0.01:1;
Uf = basismatrix(d, p, knots, w);
W = (Uf'*Uf)/101;
R = chol(W);
CM = (R*(C'*C)*R')/(h-1);
% Calculate the eigenvalues and eigenvectors of the CM
[V,D]=eig(CM);
% % Arrange feature values in descending order in DS
Ev = diag(D);
[Ev,I] = sort(Ev, 'descend');
V = V(:,I);
DS(:,1) = Ev;


% Calculate of the contribution rate
for i=1:p
    DS(i,2)=DS(i,1)/sum(DS(:,1));% single contribution rate
    DS(i,3)=sum(DS(1:i,1))/sum(DS(:,1));% cumulative contribution rate
end
% Information retention of principal components
T=0.9;
for k=1:p
    if DS(k,3) >= T
        com_num=k;
        break;
    end
end
% Extracting the eigenvectors of principal components
for j=1:com_num
    PV(:,j)=V(:,j);
end
B = R\PV;
UB = Uf*B;
UC = Uf*C';
new_score = (UB'*UC)/101;

SE = [];
YP = [];
for i = 1:1:20
    Sq = 1:1:20;
    Sq(i) = [];
    Sq_test = i;
    Sq_train = Sq;
    % parameter estimation
    X_train = X(Sq_train,:);
    [n,q] = size(X_train);
    score_train = new_score(:,Sq_train);
    score_train = score_train(:);
    theta0 = 1e-6*ones(1,q); delta0 = 1e-1; % initialization
    iterations = 100;
    
    [model] = optimization(X_train, score_train, @corrGaussian, theta0, delta0, iterations);
    
    theta0 = model.theta; delta0 = model.delta; beta0 = model.beta;
    % predict score
    ly = length(score_train);
    m = ly/n; 
    Fx = [ones(n,1),X_train];
    F = kron(Fx,ones(m,1));
    
    X_test = Xt(Sq_test,:); % using testing inputs
    Fp = [ones(1,1),X_test];
    Fp = kron(Fp,ones(m,1));
    
    Sig = eye(m); 
    Rx = corrGaussian(theta0(:), X_train);
    Sigy = kron(Rx,Sig)+delta0*eye(ly);
    Sigy_inv = eye(ly)/Sigy;
    
    R1 = corrGaussian(theta0(:), X_test, X_train);
    Cz = kron(R1,Sig);
    
    Escore = Fp*beta0+Cz*Sigy_inv*(score_train-F*beta0); % the conditional expectation of score
    
    % functional prediction
    ww = Mt(:,Sq_test);
    Ut = basismatrix(d, p, knots, ww');
    L1 = Ut*B*Escore;
    Yp = L1'.*STDt+MEANt;
    YP = [YP,Yp'];
    se = (Yp - At(Sq_test,:)).^2;
    SE = [SE;se];
end
end
