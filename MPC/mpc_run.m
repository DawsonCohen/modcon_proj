bench
exogenous

%% Luenberger Observer
L = LOB(A,C);
Ahat = [
    Ahat, zeros(n);
    -L*C, A+L*C
];
Bhat = [
    Bhat;
    zeros(n,m), -L*D, -L, B;
];
Chat = [
    zeros(4*m,n) zeros(4*m,n);
    zeros(n) eye(n);
];
Dhat = [
    Dhat;
    zeros(n,4*m);
];

%% Model Discritization
[Ad,Bd,Cd,Dd]=c2dm(Ahat,Bhat,Chat,Dhat,ts,'zoh');

n = size(Ad,1);
m = size(Bd,2);
p = size(Cd,1);

%% MPC Variable Setup
Np = 3; % prediction horizon

% Augmented state-space
Aa = [
    Ad zeros(n,p);
    Cd*Ad eye(p)
];
Ba = [Bd;Cd*Bd];
Ca = [zeros(n,p)' eye(p)];

% Outputs and reference
W = [Ca*Aa; Ca*(Aa^2); Ca*(Aa^3)];
Z = [
    Ca*Ba zeros(p,m) zeros(p,m);
    Ca*Aa*Ba Ca*Ba zeros(p,m);
    Ca*(Aa^2)*Ba Ca*Aa*Ba Ca*Ba
];

%% MPC With Constraints
N = 2000; % Time horizon
R = 0.01*diag(repmat([zeros(1,m1), ones(1,m2)],1,Np));
Q = diag(repmat([zeros(1,p1+p2), 0, 1, 0, 0, 1, 0],1,Np));
r = zeros(N+1,p);
r(:,10) = xd(:,1);
r(:,13) = xd(:,2);
del_xd = diff(xd);
del_n_proc = diff(n_proc);
del_n_sensor = diff(n_sensor);

alpha = 0.01;
beta = 0.01;

J = zeros(N+1,1);
x = zeros(N+1,n);
y = zeros(N+1,p);
u = zeros(N+1,m);
Xa = zeros(n+p,1);
nc=12; % Constraint count
del_u = zeros(Np*m,1);
G = [
    -eye(2)  zeros(2)  zeros(2);
    -eye(2)   -eye(2)  zeros(2);
    -eye(2)   -eye(2)   -eye(2);
     eye(2)  zeros(2)  zeros(2);
     eye(2)    eye(2)  zeros(2);
     eye(2)    eye(2)   -eye(2);
];
G = [zeros(nc,m1*Np) G];
mu = zeros(nc,p);

bound = [9.98; 9.98];

iters=5e3;
for k = 1:(N-Np)
    rp = r(k:k+Np-1,:);
    rp = reshape(rp',1,[])';

    for l = 1:iters % Gradient Descent
        H=[
            bound+u(k,end-1:end)';
            bound+u(k,end-1:end)'; 
            bound+u(k,end-1:end)';
            bound-u(k,end-1:end)';
            bound-u(k,end-1:end)';
            bound-u(k,end-1:end)'];
        g = G*del_u - H; % Constraint function
        gradJ = -Z'*Q*(rp - W*Xa -Z*del_u) + R*del_u;
        D_g = G;
        del_u_diff = del_u - alpha*(gradJ + D_g' * mu);
        for i = 1:Np
            del_u((i*m-m2+1):i*m) = del_u_diff((i*m-m2+1):i*m);
        end
        mu = max(mu + beta*(g),0);
        mu = max(zeros(nc,1), mu);
    end

    del_u_full = [del_xd(k,:) del_n_proc(k,:) del_n_sensor(k,:) zeros(1,m2)] + del_u(1:m)';

    J_k = 0.5*(rp-W*Xa-Z*del_u)'*Q*(rp-W*Xa-Z*del_u) + (0.5)*del_u'*R*del_u;
    u_k = u(k,:) + del_u_full;
    x_k = Ad*x(k,:)' + Bd*del_u_full';
    y_k = Ca*Xa;
    Xa = Aa*Xa+Ba*del_u_full';
    J(k+1,:) = J_k';
    x(k+1,:) = x_k';
    y(k+1,:) = y_k';
    u(k+1,:) = u_k';
end

figure(1)
plot(y(:,10), y(:,13))
hold on
plot(xd(:,1),xd(:,2))
legend(["Trajectory", "Reference"])
hold off
saveas(gcf,'mpc_lissajous.png')

fprintf("Max input: %f\n",max(u(5:6),[],'all'))
fprintf("Mean error: %f\n",mean(abs([y(:,10) y(:,13)]-xd),'all'))