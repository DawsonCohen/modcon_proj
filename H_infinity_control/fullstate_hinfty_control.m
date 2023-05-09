function F = fullstate_hinfty_control(A,B1,B2,C1,D11,D12)
    n = size(A,1);
    m1 = size(B1,2);
    p1 = size(C1,1);

    cvx_begin sdp
        variable Y(n,n) semidefinite
        variable Z(1,n)
        variable g(1)
        minimize g
        subject to
            [
                (Y*A'+A*Y+Z'*B2'+B2*Z) B1 (Y*C1'+Z'*D12') zeros(n);
                B1' -g*eye(m1) D11' zeros(m1,n);
                (C1*Y+D12*Z) D11 -g*eye(p1) zeros(p1,n)
                zeros(n) zeros(n,m1) zeros(n,1) -Y
            ] <= 0;
    cvx_end

    F = Z/Y;
end