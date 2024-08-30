function SE_RY = rmse_fixed(X,Yy,Mm,Yy1,Mm1,D,Sq_train,Sq_test)
X_train = X(Sq_train,:);
X_test = X(Sq_test,:);
[n,k] = size(X_train);
n1 = size(X_test,1);

p = 20; % the number of control points
m = 20; % dimensional of outputs
d = 3; % the order of B-spline
s = 3; % the number of replication measurements
A = (1:1:p); % the index of control points
sq = 1:60;

% Reorganize the simulation training data Y
YY = Yy(sq,Sq_train);
Y = [];
for i = 0:s-1
    YT = YY(i*m+1:i*m+m,:);
    YT = YT(:);
    Y = [Y;YT];
end

% Generate the corresponding basis matrix B
M = Mm(sq,Sq_train); MT = M(1:m,1);
knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';
Bb = basismatrix(d, p, knots, MT');
V = kron(ones(s,1),eye(n));
B = kron(V,Bb);

% Generate the simulation testing data Y1 and the corresponding basis matrix B1
Y1 = Yy1(:,Sq_test); Y1 = Y1(:);
M1 = Mm1(:,Sq_test);
B1 = [];
I = eye(n1);
for j = 0:0
    MT = M1; % M1(j*20+1:j*20+20,:)
    Bm = [];
    for i = 1:n1
        ss = MT(:,i);
        knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';
        U = basismatrix(d, p, knots, ss');
        U = kron(I(:,i),U);
        Bm = [Bm,U]; 
    end
    B1 = [B1;Bm];
end

% parameter estimation
% PAR = xlsread('result\initial_parameters.xlsx'); % initialization
% para = PAR(Sq_test,:);
% theta0 = para(1,1:k); phi0 = para(1,k+1); delta0 = para(1,k+2);
theta0 = 1e-7*ones(1,k); delta0 = 1e-3;

iterations = 100;

[model] = optimization1(X_train, X_test, Y, Y1, A, B, B1, D, @corrGaussian, s, theta0, delta0, iterations);

theta0 = model.theta; delta0 = model.delta; beta0 = model.beta;

% prediction procedure
ly = length(Y);
m = ly/(n*s); 
Fx = [ones(n,1),X_train];
f = kron(Fx,ones(m,1));
F = kron(ones(s,1),f);

Xt = X_test; % using testing inputs
Fp = [ones(n1,1),Xt];
Fp = kron(ones(51,1),Fp);

Sig = eye(p); 
Rx = corrGaussian(theta0(:), X_train);
R1 = corrGaussian(theta0(:), Xt, X_train);

Cz = kron(R1*V',Sig*Bb');
Rv = V*Rx*V'; Sigb = Bb*Sig*Bb';
[Ur,Dr] = eig(Rv); [Us,Ds] = eig(Sigb);
Urs = kron(Ur,Us); 
Drs = kron(Dr,Ds)+delta0*eye(ly); Drs = sparse(Drs);
Frs = Urs'*F; Yrs = Urs'*Y;
Crs = Cz*Urs;

Ez = Crs/Drs*(Yrs-Frs*beta0); % the conditional expectation of z
Yp = Fp*beta0+B1*Ez; % functional prediction
Yp = reshape(Yp,n1,[]);
Yt = reshape(Y1,n1,[]);
SE_RY = (Yt-Yp).^2;

end