function F = ud_8_pts(m2,m1)

m1 = [m1; ones(1, size(m1,2))];
m2 = [m2; ones(1, size(m2,2))];

L =  [];
for i=1:size(m1,2)
    L = [L; kron(m1(:,i)',m2(:,i)') ];
end

[~,~,V] = svd(L);
F = reshape(V(:,end),3,3);

%   [U,D,V] = svd(F)
%   D(end,end) = 0;
%   F = U * D * V',