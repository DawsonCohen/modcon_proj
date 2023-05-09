function norm = hinftynorm(sys)
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    
    I1 = eye(size(B',1));
    I2 = eye(size(D));

    n = size(A,1);

    eps = 1e-15;
    cvx_begin sdp quiet
        variable X(n,n) semidefinite
        variable g(1)
        minimize g
        subject to
            [
                (A'*X+X*A) X*B C';
                B'*X -g*I1 D';
                C D -g*I2
            ] <= 0;
            X >= 0;
    cvx_end

    norm = g;
end