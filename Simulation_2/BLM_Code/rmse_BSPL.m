function SE = rmse_BSPL(X,A,M,A1,M1)

h = size(X,1);
p = 15; % the number of control points
d = 2; % the order of B-spline
ss = M1(:,1);
knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';
Ut = basismatrix(d, p, knots, ss');

new_score = [];
for i = 1:h
    m1 = M(:,i);
    U = basismatrix(d, p, knots, m1');
    Y = A(:,i);
    P = U\Y;
    new_score = [new_score,P];
end

SE = [];
for i = 1:1:20
    Sq = 1:1:20;
    Sq(i) = [];
    Sq_test = i;
    Sq_train = Sq;
    % parameter estimation
    X_train = X(Sq_train,:);
    [n,p] = size(X_train);
    score_train = new_score(:,Sq_train);
    score_train = score_train(:);
    theta0 = 1e-5*ones(1,p); delta0 = 5e-3; % initialization
    iterations = 100;
    
    [model] = optimization(X_train, score_train, @corrGaussian, theta0, delta0, iterations);
    
    theta0 = model.theta; delta0 = model.delta; beta0 = model.beta;
    % predict score
    ly = length(score_train);
    m = ly/n; 
    Fx = [ones(n,1),X_train];
    F = kron(Fx,ones(m,1));
    
    X_test = X(Sq_test,:); % using testing inputs
    Fp = [ones(1,1),X_test];
    Fp = kron(Fp,ones(m,1));
    %
    Sig = eye(m); 
    Rx = corrGaussian(theta0(:), X_train);
    Sigy = kron(Rx,Sig)+delta0*eye(ly);
    Sigy_inv = eye(ly)/Sigy;
    
    R1 = corrGaussian(theta0(:), X_test, X_train);
    Cz = kron(R1,Sig);
    
    Escore = Fp*beta0+Cz*Sigy_inv*(score_train-F*beta0); % the conditional expectation of score
    
    % functional prediction
    Yp = Ut*Escore;
    se = (Yp - A1(:,Sq_test)).^2;
    SE = [SE;se'];
end
end
