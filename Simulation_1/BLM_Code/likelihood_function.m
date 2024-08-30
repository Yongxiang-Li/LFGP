function z = likelihood_function(X, Y, corr_2, theta, delta)

n = size(X,1); % the number of inputs
ly = length(Y);
m = ly/n; % m-dimensional outputs corresponding to each input
Fx = [ones(n,1),X];
F = kron(Fx,ones(m,1));
Ft = kron(Fx',ones(1,m)); % transpose of F
pq = size(F,2);

Sig = eye(m);
Rx = corr_2(theta(:), X); % the correlation matrix of inputs
Sigy = kron(Rx,Sig)+delta*eye(ly);
Sigy_inv = eye(ly)/Sigy;
triangle = chol(Sigy);
diagonal = diag(triangle);

log_Sigy_det = 2*sum(log(diagonal)); 

FK = Ft*Sigy_inv*F+1e-5*eye(pq);

beta = FK\Ft*Sigy_inv*Y;
sigma = (Y'-beta'*Ft)*Sigy_inv*(Y-F*beta)/ly;

z = (ly*log(sigma)+log_Sigy_det)/2+ly/2; % negative log likelihood

end