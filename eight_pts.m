function F = eight_pts(m2, m1)

    %m1 e m2 devono avere la stessa dimensione, lo diamo per scontato
    m1 = [m1; ones(1, size(m1, 2))]; %per rendere le coordinate omogenee  
    m2 = [m2; ones(1, size(m2, 2))];
    
    L = [];
    for i = 1: size(m1,2)
        L = [L;  kron( m1(:,i)' , m2(:,i)') ];
    end
    
    [~,~,V] = svd(L);
    F = reshape(V(:,end),3,3);
    
end
