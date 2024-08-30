function [SE,YP] = rmse_BSPL_1(X,Xt,A,At,M,Mt)

h = size(X,1);
p = 15; % % the number of control points
d = 2; % the order of B-spline
knots = [zeros(1,d) linspace(0,1,p-d+1) ones(1,d)]';

new_score = [];
for j = 1:h
    ss = M(:,j);
    Y = A(:,j);
    U = basismatrix(d, p, knots, ss');
    P = U\Y;
    new_score = [new_score,P];
end

SE = [];
YP = [];
for i = 1:1:25
    Sq = 1:1:30;
    Sq_test = i;
    Sq_train = Sq;
    % parameter estimation
    X_train = X(Sq_train,:);
    [n,q] = size(X_train);
    score_train = new_score(:,Sq_train);
    score_train = score_train(:);
    theta0 = 1e-5*ones(1,q); delta0 = 5e-2; % initialization
    iterations = 100;
    tic;
    [model] = optimization_1(X_train, score_train, @corrGaussian, theta0, delta0, iterations);
    t = toc;
    fprintf('the total running time is %7.4f s\n',t);
    theta0 = model.theta; delta0 = model.delta; beta0 = model.beta;
    % predict score
    ly = length(score_train);
    m = ly/n; 
    Fx = [ones(n,1),X_train];
    F = kron(Fx,ones(m,1));
    
    X_test = Xt(Sq_test,:); % using testing inputs
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
    ss = Mt(:,Sq_test);
    U = basismatrix(d, p, knots, ss');
    Yp = U*Escore;
    YP = [YP,Yp];
    se = (Yp - At(:,Sq_test)).^2;
    SE = [SE;se'];
end
end
