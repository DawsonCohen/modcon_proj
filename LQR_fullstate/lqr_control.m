bench

sys = ss(A,B,C,D);
sysd = c2d(sys,ts);

Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

dt = 1e-5;
tf = 1e-3;
t = 0:dt:tf;
td = 0:ts:tf;
N = length(t);
Nd = length(td);

xd_step = 1e-3;
x0 = [0 0 0 0 0 0]';
xd = [0 1 0 0 1 0]'*xd_step;
u = repmat(xd',N,1);
ud = repmat(xd',Nd,1);

eps = 100000;
delta = 1e-1;
Q = eps*diag([0 1 0 0 1 0]);
% Q = eps*eye(6);
R = delta*eye(2);

K = clqr(A,B,Q,R);
Kd = dlqr(Ad,Bd,Q,R);

Acl = A+B*K;
Bcl = -B*K;
Ccl = C;
Dcl = [];

% Acld = Ad+Bd*Kd;
% Bcld = -Bd*Kd;
% Ccld = Cd;
% Dcld = [];

syscl = ss(Acl, Bcl, Ccl, Dcl);
% syscld = ss(Acld, Bcld, Ccld, Dcld, ts);

[mag, phase, wout] = bode(syscl);

% fig = figure();
% subplot(4,2,1)
% semilogx(wout,reshape(mag(1,2,:),[],1))
% ylabel('In: V_x')
% title("Out: \theta_x",'FontWeight','Normal')
% xlim()
% subplot(4,2,2)
% semilogx(wout,reshape(mag(1,5,:),[],1))
% title("Out: \theta_y",'FontWeight','Normal')
% 
% subplot(4,2,3)
% semilogx(wout,reshape(mag(2,2,:),[],1))
% ylabel("In: V_y")
% subplot(4,2,4)
% semilogx(wout,reshape(mag(2,5,:),[],1))
% 
% subplot(4,2,5)
% semilogx(wout,reshape(phase(1,2,:),[],1))
% ylabel("In: V_x")
% subplot(4,2,6)
% semilogx(wout,reshape(phase(1,5,:),[],1))
% subplot(4,2,7)
% semilogx(wout,reshape(phase(2,2,:),[],1))
% ylabel("In: V_y")
% subplot(4,2,8)
% semilogx(wout,reshape(phase(2,5,:),[],1))
% 
% sgtitle("Bode Plot for \theta_x and \theta_y")
% 
% han=axes(fig,'visible','off'); 
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% 
% Xlb=-1;                    % set horizontally at midpoint
% Ylb=mean(ylim);            % and just 1% below minimum y value
% ylabel(han,'Phase (deg) ; Magnitude (dB)','Position', [-0.11, 0.5]);
% xlabel(han,'frequency (rad/s)');
% 
% saveas(gcf,"./figures/lqr_bode.png")
% 
% % % [x,t,u] = step(syscl);
[y,t,x] = lsim(syscl,u,t);

figure()
subplot(2,1,1)
plot(t,x(:,2)*1e3);
hold on
plot(t,u(:,2)*1e3, "--");
title("LQR Step Response")
legend("response", "target")
% xlim([0 10])
ylim([0,xd_step*1e3*1.2])
ylabel("\theta_x (mrad)")
hold off

subplot(2,1,2)
plot(t,x(:,5)*1e3);
hold on
plot(t,u(:,5)*1e3, "--");
% xlim([0 10])
ylim([0,xd_step*1e3*1.2])
xlabel("time (s)")
ylabel("\theta_y (mrad)")
hold off

% saveas(gcf,"figures/lqr_step.png")