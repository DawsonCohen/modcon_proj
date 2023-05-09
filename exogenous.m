bench

tf = 1;         % T final
t = 0:ts:tf;    % time vector
N = length(t); 

%% Lissajous Command
x = 3*sin(8*pi*(t/tf));
y = 4*sin(6*pi*(t/tf));

xd_sig = [x' y'];
xd = [x' y'] * 1e-3;

plot(xd(:,1),xd(:,2));
saveas(gcf, 'lissajous.png')

%% Process Noise
sig_proc = 1e-4;
n_proc = randn(N,m)*sig_proc;

%% Sensor Noise
sig_dfsm = 1e-1;
f_dfsm = 180;
N_dfsm = ceil((f_dfsm/fs)*N);
noise = randn(N_dfsm,m)*sig_dfsm;
n_dfsm = resample(noise, fs, f_dfsm);
n_dfsm = n_dfsm(1:N,:);

sig_signal = 1e-3;
n_signal = randn(N,m)*sig_signal;

n_sensor = n_dfsm + n_signal;

u = [xd_sig n_proc n_sensor];

x0 = [0 0 0 0 0 0]';
