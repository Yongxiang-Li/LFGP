function z = likelihood(X, Xt, Y, Y1, A, U, U1, D, corr_2, s, theta, delta)
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
q = size(D,1);


Sig = eye(p); % the correlation matrix of control points
Rx = corr_2(theta(:), X); % the correlation matrix of inputs
R1 = corrGaussian(theta(:), Xt, X);
Rt = kron(R1,Sig);
Sig_inv = eye(p)/Sig;
Rx_inv = eye(n)/Rx;
R_inv = kron(Rx_inv,Sig_inv);
omega = delta*R_inv+U'*U + 1e-5*eye(n*p);
omega_inv = eye(n*p)/omega;
% omega1 = 1e10*(delta*R_inv+U'*U) + 1e-5*eye(n*p);
% omega_inv = 1e10*eye(n*p)/omega1;
triangle = chol(omega);
diagonal = diag(triangle);

Sig_det=det(Sig);
Rx_det=det(Rx);
log_omega_det = 2*sum(log(diagonal));

Sigy_inv = eye(ly)/delta-U*omega_inv*U'/delta; % SMW
log_Sigy_det = (ly-k)*log(delta) + p*log(Rx_det) + n*log(Sig_det) + log_omega_det; % sylvester determinant
gamma = D*Rt*U'*Sigy_inv;
FK = Ft*Sigy_inv*F+1e-5*eye(pq);
Fg = gamma*F/FK*Ft;
FG = Fg*gamma';

[Q,S] = qr(FG);
S = 1e10*S+1e-5*eye(q);
K = Q'*(Fg*Sigy_inv*Y-gamma*Y);
lamda = 1e10*eye(q)/S*K;
% lamda = FG\(Fg*Sigy_inv*Y-gamma*Y);

beta = FK\Ft*(Sigy_inv*Y-gamma'*lamda);
sigma = (Y'-beta'*Ft)*Sigy_inv*(Y-F*beta)/ly;

z = (ly*log(sigma)+log_Sigy_det)/2+ly/2; % negative log likelihood

end