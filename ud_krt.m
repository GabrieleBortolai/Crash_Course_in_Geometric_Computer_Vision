function [K,R,t] = ud_krt(P)

[Q,U] = qr(inv(P(1:3,1:3)));

D = diag(sign(diag(U)) .* [-1; -1; 1]);
Q = Q*D;
U = D*U;

s = det(Q);
R = s*Q';
t = s * U * P(1:3,4);
K = inv(U./U(3,3));
