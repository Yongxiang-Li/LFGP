function z = objective_function(X, Xt, Y, Y1, A, U, U1, D, corr_1, corr_2, s, theta, phi, delta)
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


Sig = corr_1(phi, A); % the correlation matrix of control points
Rx = corr_2(theta(:), X); % the correlation matrix of inputs
R1 = corrGaussian(theta(:), Xt, X);
Rt = kron(R1,Sig);
Sig_inv = eye(p)/Sig;
Rx_inv = eye(n)/Rx;
R_inv = kron(Rx_inv,Sig_inv);
omega = 1e10*(delta*R_inv+U'*U);
triangle = chol(omega);
tri_inv = eye(n*p)/triangle;
omega_inv = 1e10*(tri_inv*tri_inv');


Sigy_inv = eye(ly)/delta-U*omega_inv*U'/delta; % SMW
gamma = D*Rt*U'*Sigy_inv;
FK = Ft*Sigy_inv*F+1e-5*eye(pq);
Fg = gamma*F/FK*Ft;
FG = Fg*gamma';

[Q,S] = qr(FG);
S = 1e10*S+1e-5*eye(q);
K = Q'*(Fg*Sigy_inv*Y-gamma*Y);
lamda = 1e10*eye(q)/S*K;
beta = FK\Ft*(Sigy_inv*Y-gamma'*lamda);

n1 = size(Xt,1);
Fp = [ones(n1,1),Xt];
Ez = Rt*U'*Sigy_inv*(Y-F*beta);
Fp = kron(ones(51,1),Fp);
Yp = Fp*beta+U1*Ez;
Yp = reshape(Yp,n1,[]);
Yt = reshape(Y1,n1,[]);
SE_RY = (Yt-Yp).^2;
z = sqrt(mean(SE_RY));
end


