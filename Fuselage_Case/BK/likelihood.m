function z = likelihood(X, Y, V, corr_2, theta, delta)

n = size(X,1); % the number of inputs
ly = length(Y);
F = V;
Ft = F'; % transpose of F
q = size(F,2);

Rx = corr_2(theta(:), X); % the correlation matrix of inputs
Rx = Rx + delta*eye(n);
% Sigy_inv = eye(n)/Rx;

triangle = chol(Rx);
tr_inv = eye(n)/triangle;
Sigy_inv = tr_inv*tr_inv';

diagonal = diag(triangle);
log_det = sum(2*log(diagonal));

FK = Ft*Sigy_inv*F+1e-5*eye(q);
beta = FK\Ft*Sigy_inv*Y;
sigma = (Y'-beta'*Ft)*Sigy_inv*(Y-F*beta)/ly;
z = (ly*log(sigma)+log_det)/2+ly/2; % negative log likelihood

end