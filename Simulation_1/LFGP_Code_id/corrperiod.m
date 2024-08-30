function [ C ] = corrperiod( phi, A )
% period correlation between two control points with parameter phi
% A is a row vector, representing the index of control points
% phi is column vector
n = length(A);
K = A'*ones(1,n);
Phi = pi*(K-K')/n;
Phis = sin(Phi).^2;
C = exp(-phi*Phis);
C = C + 1e-8*eye(n);
end