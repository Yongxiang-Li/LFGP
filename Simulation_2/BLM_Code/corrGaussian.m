function [ C ] = corrGaussian( theta, X, Y )
% Gaussian correlation between two points x and y with parameter theta
% x and y are row vectors
% theta is column vector
    n = size(X,1);
    X = X.*(ones(n,1)*sqrt(theta'));
    if nargin == 2
        XY = X*X';
        C = exp(2*XY-(diag(XY)*ones(1,n)+ones(n,1)*diag(XY)'));
        C = C + 1e-8*eye(size(C));
    else
        k = size(Y,1);
        Y = (Y.*(ones(k,1)*sqrt(theta')))';
        C = exp(2*X*Y-((X.^2)*ones(size(theta))*ones(1,k)+ones(n,1)*ones(size(theta'))*(Y.^2)));
    end
end