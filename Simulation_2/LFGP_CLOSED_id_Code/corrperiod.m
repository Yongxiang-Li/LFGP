function [ C ] = corrperiod( phi, A, B )
% period correlation between two control points with parameter phi
% A is a row vector, representing the index of control points
% phi is column vector
n = length(A);
if nargin == 2
    K = A'*ones(1,n);
    Phi = pi*(K-K')/(n-1);
    Phis = sin(Phi).^2;
    C = exp(-phi*Phis);
    C = C + 1e-8*eye(n);
else
    n1 = length(B);
    K = A'*ones(1,n1);
    K1 = ones(n,1)*B;
    Phi = pi*(K-K1)/360;
    Phis = sin(Phi).^2;
    C = exp(-phi*Phis);
end
end