function SE_RY = rmse_fixed(X,Yy,Mm,Yy1,Mm1,Sq_train,Sq_test)
X_train = X(Sq_train,:);
X_test = X(Sq_test,:);
n = size(X_train,1);
n1 = size(X_test,1);

p = 20; % the number of control points
m = 20; % dimensional of outputs
d = 3; % the order of B-spline
s = 3; % the number of replication measurements
A = (1:1:p); % the index of control points
sq = 1:1:60;

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
B = basismatrix(d, p, knots, MT');

% Generate the simulation testing data Y1 and the corresponding basis matrix B1
Y1 = Yy1(:,Sq_test); Y1 = Y1(:);
M1 = Mm1(:,Sq_test);
B1 = [];
I = eye(n1);
for j = 0:0
    MT = M1; %M1(j*m+1:j*m+m,:);
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
theta0 = [1e-7,1e-3]; phi0 = 4e4; delta0 = 1e-3; % initialization
iterations = 100;

[model] = optimization_fix(X_train, Y, A, B, @correxponent, @corrGaussian, s, theta0, phi0, delta0, iterations);

theta0 = model.theta; phi0 = model.phi; delta0 = model.delta; beta0 = model.beta; % sigma0 = model.sigma;

% prediction procedure
ly = length(Y);
m = ly/(n*s); % m-dimensional outputs corresponding to each input
Fx = [ones(n,1),X_train];
f = kron(Fx,ones(m,1));
F = kron(ones(s,1),f);

Xt = X_test; % using testing inputs
Fp = [ones(n1,1),Xt];

Sig = correxponent(phi0, A); 
Rx = corrGaussian(theta0(:), X_train);
R1 = corrGaussian(theta0(:), Xt, X_train);
V = kron(ones(s,1),eye(n));

Cz = kron(R1*V',Sig*B');
Rv = V*Rx*V'; Sigb = B*Sig*B';
[Ur,Dr] = eig(Rv); [Us,Ds] = eig(Sigb);
Urs = kron(Ur,Us); 
Drs = kron(Dr,Ds)+delta0*eye(ly); Drs = sparse(Drs);
Frs = Urs'*F; Yrs = Urs'*Y;
Crs = Cz*Urs;

Ez = Crs/Drs*(Yrs-Frs*beta0); % the conditional expectation of z
Fp = kron(ones(51,1),Fp);
Yp = Fp*beta0+B1*Ez; % functional prediction
Yp = reshape(Yp,n1,[]);
Yt = reshape(Y1,n1,[]);
SE_RY = (Yt-Yp).^2;
end