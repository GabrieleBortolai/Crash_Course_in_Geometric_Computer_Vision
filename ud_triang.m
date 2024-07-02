function M = ud_triang(P,m)

L = [];

for j=1:length(P)
   L = [L; skew([m{j};1])*P{j}];
end

[U,D,V] = svd(L);

M = V(:,end);

M = M./M(end);
M(end) = [];


