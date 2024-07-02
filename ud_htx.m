function Y = ud_htx(T, X)

Y = T *  [X; ones(1,size(X,2))];
nr = size(Y,1);
Y = Y./repmat( Y(nr,:), nr, 1);
Y(nr,:)= [];
