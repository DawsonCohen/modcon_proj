bench
% load bench_ABCD.mat

sys = ss(A,B,C,D);
dt = 0.0001;

% [y,t,x] = lsim(sys,u,t,x0);
t = 0:dt:1e-3;
N = length(t);

x0 = [0 0 0 0 0 0]';
x = zeros(6,length(T));
x(:,1) = x0;
ui = [10 10]';
u = repmat(ui',N,1);

figure()
step(sys)
saveas(gcf,"./figures/open_step.png")
figure()
impulse(sys)
saveas(gcf,"./figures/open_impulse.png")

% [y,t,x] = lsim(sys,u,t)

% figure(1)
% plot(t,x(:,2));
% figure(2)
% plot(t,x(:,5));