function [SE_RY,Yp,model] = rmse_random(X,Xxt,Yy,Yyt,Mm,Mmt,D,Sq_train,Sq_test)
X_train = X(Sq_train,:);
X_test = Xxt(Sq_test,:);
[n,k] = size(X_train);
n1 = size(X_test,1);

p = 20; % the number of control points 15
m = 20; % dimensional of outputs
d = 3; % the order of B-spline 5
s = 3; % the number of replication measurements
A = (1:1:p); % the index of control points

% Reorganize the simulation training data Y
YY = Yy(:,Sq_train);
Y = [];
for i = 0:s-1
    YT = YY(i*m+1:i*m+m,:);
    YT = YT(:);
    Y = [Y;YT];
end

% Generate the simulation training data Y and the corresponding basis matrix B
M = Mm(:,Sq_train);
B = [];
I = eye(n);
for j = 0:s-1
    MT = M(j*m+1:j*m+m,:);
    Bm = [];
    for i = 1:n
        ss = MT(:,i);
        knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';
        U = basismatrix(d, p, knots, ss');
        U = kron(I(:,i),U);
        Bm = [Bm,U]; 
    end
    B = [B;Bm];
end
% Generate the simulation testing data Y1 and the corresponding basis matrix B1
Y1 = Yyt(:,Sq_test); Y1 = Y1(:);
M1 = Mmt(:,Sq_test);
B1 = [];
I = eye(n1);
for j = 0:s-1
    MT = M1(j*m+1:j*m+m,:);
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
theta0 = 5e-7*ones(1,k); phi0 = 200; delta0 = 1e-3;

iterations = 100;

[model] = optimization(X_train, X_test, Y, A, B, D, @corrperiod, @corrGaussian, s, theta0, phi0, delta0, iterations);

theta0 = model.theta; phi0 = model.phi; delta0 = model.delta; beta0 = model.beta;

% prediction procedure
ly = length(Y);
m = ly/(n*s); 
Fx = [ones(n,1),X_train];
f = kron(Fx,ones(m,1));
F = kron(ones(s,1),f);

Xt = X_test; % using testing inputs
Fp = [ones(n1,1),Xt];
Fp = kron(ones(60,1),Fp);

Sig = corrperiod(phi0, A); 
Rx = corrGaussian(theta0(:), X_train);
R1 = corrGaussian(theta0(:), Xt, X_train);
Sig_inv = eye(p)/Sig;
Rx_inv = eye(n)/Rx;
R_inv = kron(Rx_inv,Sig_inv);

Cz = kron(R1,Sig)*B';
omega = delta0*R_inv+B'*B;
triangle = chol(omega);
Omega = triangle'*triangle;
omega_inv = eye(n*p)/Omega;
L = chol(omega_inv)';

Yb = L'*B'*Y; Fb = L'*B'*F;
Ey = Y-F*beta0; Eb = Yb-Fb*beta0;
Cb = Cz*B*L;

Ez = (Cz*Ey-Cb*Eb)/delta0; % the conditional expectation of z
Yp = Fp*beta0+B1*Ez; % functional prediction
Yp = reshape(Yp,n1,[]);
Yt = reshape(Y1,n1,[]);
SE_RY = (Yt-Yp).^2;

end