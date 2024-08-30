function [model] = optimization_1(X, Y, corr_2, theta0, delta0, iterations)
% Y is smn functional outputs for s replications
Y = Y(:);
[n,q] = size(X);
ly = length(Y);
m = ly/n; 
Fx = [ones(n,1),X];
F = kron(Fx,ones(m,1));
Ft = kron(Fx',ones(1,m)); 
pq = size(F,2);

maxLS0 = -likelihood_function(X, Y, corr_2, theta0, delta0);
dLS=abs(maxLS0);
options=optimset('Algorithm','sqp', 'Display','off');
i=0;
while dLS>1e-3 && i<iterations    %(i can be smaller otherwise the convergence is time-consuming)
    i=i+1; 
    % MLE to parameter estimation
    theta0=fmincon(@(theta) likelihood_function(X, Y, corr_2, theta, delta0), ...
        theta0, [], [], [], [], 1e-8*ones(1,q), 6e-3*ones(1,q), [], options); 
    delta0=fmincon(@(delta) likelihood_function(X, Y, corr_2, theta0, delta), ...
        delta0, [], [], [], [], 1e-4, 1e-1, [], options); 
    % calculate likelihood
    maxLS = -likelihood_function(X, Y, corr_2, theta0, delta0);
    dLS=maxLS-maxLS0; 
    maxLS0=maxLS;
end

Sig = eye(m);
Rx = corr_2(theta0(:), X); % the correlation matrix of inputs
Sigy = kron(Rx,Sig)+delta0*eye(ly);
Sigy_inv = eye(ly)/Sigy;
FK = Ft*Sigy_inv*F+1e-5*eye(pq);

beta0 = FK\Ft*Sigy_inv*Y;
sigma0 = (Y'-beta0'*Ft)*Sigy_inv*(Y-F*beta0)/ly;

model.theta = theta0;
model.delta = delta0;
model.beta = beta0;
model.sigma = sigma0;
model.iteration = i;

end