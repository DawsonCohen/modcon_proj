bench

P = lyap(A,B*B');
min(eig(P))

W = lyap(A',C'*C);
min(eig(W))