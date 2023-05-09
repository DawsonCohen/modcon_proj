function L = LOB(A,C)
    n = size(A,1);
    p = size(C,1);
    
    eps = 1e-9
    cvx_begin sdp
        variable W(n,n) symmetric
        variable V(p,n)
        A'*W + W*A + C'*V + V'*C <= -eps*eye(n)
        W >= eps*eye(n)
    cvx_end
    
    L = W\(V');
end