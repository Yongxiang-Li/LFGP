function z = likelihood_random(X, Y, A, U, corr_1, corr_2, s, theta, phi, delta)
% U is a smn*nk matrix of basis functions, m denotes the hierarchy of output Y and k denotes the hierarchy of control points
n = size(X,1); % the number of inputs
ly = length(Y);
p = length(A);
k = size(U,2);
m = ly/(n*s); % m-dimensional outputs corresponding to each input
Fx = [ones(n,1),X];
f = kron(Fx,ones(m,1));
ft = kron(Fx',ones(1,m)); % transpose of f
F = kron(ones(s,1),f);
Ft = kron(ones(1,s),ft); % transpose of F
pq = size(F,2);

Sig = corr_1(phi, A); % the correlation matrix of control points
Rx = corr_2(theta(:), X); % the correlation matrix of inputs
Sig_inv = eye(p)/Sig;
Rx_inv = eye(n)/Rx;
R_inv = kron(Rx_inv,Sig_inv);
omega = delta*R_inv+U'*U;
triangle = chol(omega);
Omega = triangle'*triangle;
omega_inv = eye(n*p)/Omega;

L = chol(omega_inv)';
diagonal = diag(triangle);

Sig_det=det(Sig);
Rx_det=det(Rx);
log_omega_det=2*sum(log(diagonal));
log_Sigy_det = (ly-k)*log(delta) + p*log(Rx_det) + n*log(Sig_det) + log_omega_det; % sylvester determinant

Yb = L'*U'*Y; Fb = L'*U'*F; Fbt = Ft*U*L;
FK = Ft*F-Fbt*Fb+1e-5*eye(pq);
beta = FK\(Ft*Y-Fbt*Yb);
sigma = (norm(Y-F*beta)^2-norm(Yb-Fb*beta)^2)/(delta*ly);

z = (ly*log(sigma)+log_Sigy_det)/2+ly/2; % negative log likelihood
end


