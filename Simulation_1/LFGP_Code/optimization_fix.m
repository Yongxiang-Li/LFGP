function [model] = optimization_fix(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0, iterations)
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

maxLS0 = -likelihood_fix(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0);
dLS=abs(maxLS0);
options=optimset('Algorithm','sqp','Display','off','MaxFunEvals',1000,'MaxIter',50);
i=0;
while dLS>1e-3 && i<iterations    %(i can be smaller otherwise the convergence is time-consuming)
    i=i+1; 
    % MLE to parameter estimation
    theta0=fmincon(@(theta) likelihood_fix(X, Y, A, U, corr_1, corr_2, s, theta, phi0, delta0), ...
        theta0, [], [], [], [], [1e-8,1e-8], [1e-2,1e-2], [], options);
    phi0=fmincon(@(phi) likelihood_fix(X, Y, A, U, corr_1, corr_2, s, theta0, phi, delta0), ...
        phi0, [], [], [], [], 4e3, 1.2e5, [], options);
    delta0=fmincon(@(delta) likelihood_fix(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta), ...
        delta0, [], [], [], [], 1e-4, 1e-1, [], options);
    % calculate likelihood
    maxLS = -likelihood_fix(X, Y, A, U, corr_1, corr_2, s, theta0, phi0, delta0);
    dLS=maxLS-maxLS0; 
    maxLS0=maxLS;
end

Sig = corr_1(phi0, A); % the correlation matrix of control points
Rx = corr_2(theta0(:), X); % the correlation matrix of inputs
V = kron(ones(s,1),eye(n));

Rv = V*Rx*V'; Sigb = U*Sig*U';
[Ur,Dr] = eig(Rv); [Us,Ds] = eig(Sigb);
Urs = kron(Ur,Us); 
Drs = kron(Dr,Ds)+delta0*eye(ly); Drs = sparse(Drs);
Frs = Urs'*F; Frst = Ft*Urs; Yrs = Urs'*Y;
FK = Frst/Drs*Frs+1e-5*eye(pq);

beta0 = FK\Frst/Drs*Yrs;
% sigma0 = (Yrs'-beta0'*Frst)/Drs*(Yrs-Frs*beta0)/ly;

model.theta = theta0;
model.phi = phi0;
model.delta = delta0;
model.beta = beta0;
% model.sigma = sigma0;
model.iteration = i;

end