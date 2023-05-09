function [Ak, Bk, Ck, Dk, Ak2, Bk2, Ck2, Dk2] = hinfty_control(A,B1,B2,C1,C2,D11,D12,D21,D22)
    n = size(A,1);
    m1 = size(B1,2);
    m2 = size(B2,2);
    p1 = size(C1,1);
    p2 = size(C2,1);

    eps = 1e-15;
    cvx_begin sdp
        variable X1(n,n) semidefinite
        variable Y1(n,n) semidefinite
        variable An(n,n)
        variable Bn(n,p2)
        variable Cn(m2,n)
        variable Dn(m2,p2)
        variable g(1)
        minimize g
        subject to
            [
                Y1 eye(n);
                eye(n) X1;
            ] >= eps*eye(2*n);
            [
                (A*Y1 + Y1*A' + B2*Cn + Cn'*B2'), (A' + An + (B2*Dn*C2)')', (B1+B2*Dn*D21), (C1*Y1+D12*Cn)';
                (A' + An + (B2*Dn*C2)'), ((X1*A) + A'*X1 + Bn*C2 + C2'*Bn'), (X1*B1+Bn*D21),  (C1+D12*Dn*C2)';
                (B1+B2*Dn*D21)', (X1*B1+Bn*D21)', -g*eye(m1), (D11+D12*Dn*D21)';
                (C1*Y1+D12*Cn), (C1+D12*Dn*C2), (D11+D12*Dn*D21), -g*eye(p1);
            ] <= -eps*eye(2*n+m1+p1);
    cvx_end

    X2 = eye(n)-X1*Y1;
    Y2 = eye(n);

    cvx_begin quiet
        variable Ak2(n,n)
        variable Bk2(n,p2)
        variable Ck2(p2,n)
        variable Dk2(m2,p2)
        
        [
            Ak2 Bk2;
            Ck2 Dk2
        ] == [
            X2 X1*B2;
            zeros(m2,n) eye(p2);
        ]\( ...
            [
                An Bn;
                Cn Dn;
            ] - [
                X1*A*Y1 zeros(n,p2);
                zeros(m2,n) zeros(m2,p2);
            ]) / ...
        [
            Y2' zeros(n,p2);
            C2*Y1 eye(p2);
        ]
    cvx_end

    Dk = (eye(m2) + Dk2*D22)\Dk2;
    Bk = Bk2 * (eye(p2) - D22*Dk);
    Ck = (eye(m2)-Dk*D22)*Ck2;
    Ak = Ak2 - Bk/(eye(p2)-D22*Dk)*D22*Ck;
end