function [model] = optimization_random(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0, iterations)
% Y is smn functional outputs for s replications and U is smn*nk basis matrix
Y = Y(:);
n = size(X,1);
p = length(A);
ly = length(Y);
m = ly/(n*s); % m-dimensional outputs corresponding to each input
Fx = [ones(n,1),X];
f = kron(Fx,ones(m,1));
ft = kron(Fx',ones(1,m)); % transpose of f
F = kron(ones(s,1),f);
Ft = kron(ones(1,s),ft); % transpose of F
pq = size(F,2);

maxLS0 = -likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0);
dLS=abs(maxLS0);
options=optimset('Algorithm','sqp', 'Display','off');
i=0;
while dLS>1e-3 && i<iterations    %(i can be smaller otherwise the convergence is time-consuming)
    i=i+1; 
    % MLE to parameter estimation
    theta0=fmincon(@(theta) likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta, phi0, delta0), ...
        theta0, [], [], [], [], [1e-8,1e-8], [1e-2,1e-2], [], options);
    phi0=fmincon(@(phi) likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta0, phi, delta0), ...
        phi0, [], [], [], [], 4e3, 1.2e5, [], options);
    delta0=fmincon(@(delta) likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta), ...
        delta0, [], [], [], [], 1e-4, 1e-1, [], options);
    % calculate likelihood
    maxLS = -likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0);
    dLS=maxLS-maxLS0; 
    maxLS0=maxLS;
end

Sig = corr_1(phi0, A); % the correlation matrix of control points
Rx = corr_2(theta0(:), X); % the correlation matrix of inputs
Sig_inv = eye(p)/Sig;
Rx_inv = eye(n)/Rx;
R_inv = kron(Rx_inv,Sig_inv);
omega = delta0*R_inv+U'*U;
triangle = chol(omega);
Omega = triangle'*triangle;
omega_inv = eye(n*p)/Omega;
L = chol(omega_inv)';

Yb = L'*U'*Y; Fb = L'*U'*F; Fbt = Ft*U*L;
FK = Ft*F-Fbt*Fb+1e-5*eye(pq);

beta0 = FK\(Ft*Y-Fbt*Yb);
% sigma0 = (norm(Y-F*beta0)^2-norm(Yb-Fb*beta0)^2)/(delta0*ly);

model.theta = theta0;
model.phi = phi0;
model.delta = delta0;
model.beta = beta0;
% model.sigma = sigma0;
model.iteration = i;

end