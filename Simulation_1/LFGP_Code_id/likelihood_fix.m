function z = likelihood_fix(X, Y, A, U, corr_2, s, theta, delta)
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

Sig = eye(p); % the correlation matrix of control points is identity matrix
Rx = corr_2(theta(:), X); % the correlation matrix of inputs
V = kron(ones(s,1),eye(n));

Rv = V*Rx*V'; Sigb = U*Sig*U';
[Ur,Dr] = eig(Rv); [Us,Ds] = eig(Sigb);
Urs = kron(Ur,Us); 
Drs = kron(Dr,Ds)+delta*eye(ly); Drs = sparse(Drs);
diagD = diag(Drs);
Frs = Urs'*F; Frst = Ft*Urs; Yrs = Urs'*Y;
FK = Frst/Drs*Frs+1e-5*eye(pq);

beta = FK\Frst/Drs*Yrs;
sigma = (Yrs'-beta'*Frst)/Drs*(Yrs-Frs*beta)/ly;
z = (ly*log(sigma)+sum(log(diagD)))/2+ly/2; % negative log likelihood

end
