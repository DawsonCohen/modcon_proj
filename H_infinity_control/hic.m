bench
exogenous

A = Ahat;
B1 = B1hat;
B2 = B2hat;
C1 = C1hat;
C2 = C2hat;
D11 = D11hat;
D12 = D12hat;
D21 = D21hat;
D22 = D22hat;

[Ak, Bk, Ck, Dk] = hinfty_control( ...
    A,B1,B2,C1,C2,D11,D12,D21,D22);

Q = eye(p2)/(eye(p2)-D22*Dk);
Acl = [
    A+B2*Dk*Q*C2, B2*(eye(m2)+Dk*Q*D22)*Ck;
    Bk*Q*C2, Ak+Bk*Q*D22*Ck;
];
Bcl = [
    B1+B2*Dk*Q*D21;
    Bk*Q*D21;
];
Ccl = [
    C1+D12*Dk*Q*C2, D12*(eye(m2)+Dk*Q*D22)*Ck
];
Dcl = D11+D12*Dk*Q*D21;

syscl = ss(Acl, Bcl, Ccl, Dcl);

[y,t,x] = lsim(syscl,u,t);
figure(3)
plot(y(:,1),y(:,2));
hold on;
plot(xd_sig(:,1),xd_sig(:,2))
legend(["Trajectory", "Reference"])
hold off
saveas(gcf,'hic_liss.png')
figure(4)
lsim(syscl,u,t);
saveas(gcf,'hic_fullstate.png')

err = mean(abs(y(:,1:2)-xd_sig));
figure(5)
sigma(syscl)

saveas(gcf,'hic_sigma.png')