function [s, R, t] = ud_opa(D,M)

n = size(M,2);
J = eye(n) - ones(n,1) * ones(1,n)/n;

[U,~,V] = svd(D*J*M');

R = U*diag([1,1,det(V'*U)])*V';
s = trace(D*J*M'*R')/trace(M*J*M');
t = (D/s - R*M)*ones(n,1)/n;


