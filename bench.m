step_size = 0.001;

T = 67e-6;      % FSM response time
Dm = 1.245;     % Distance to Mirror - m
zeta_x = 0.9;   % FSM Damping Ratio, x
zeta_y = 0.9;   % FSM Damping Ratio, y
omega_x = 5655; % FSM natural frequency, x
omega_y = 5184; % FSM natural frequency, y

Gfx = 2000;     % FPA callibration, x - px / rad
Gfy = 2000;     % FPA callibration, y - px / rad
Gmx = 52.4e-3;  % Mirror callibration, x - rad / V
Gmy = 52.4e-3;  % Mirror callibration, x - rad / V
Ax = -2e-2;     % Cross-coupling factor, x
Ay = -2e-2;     % Cross-coupling factor, y

Vsupply = 10;   % Maximum input

fs = 2e3;       % Sample rate
ts = 1/fs;      % Sample period

A = [
    -1/T Gfx*Dm/T 0 0 0 0;
    0 0 1 0 0 0;
    0 -omega_x^2 -2*zeta_x*omega_x 0 Ax*omega_x^2 Ax*2*zeta_x*omega_x;
    0 0 0 -1/T Gfy*Dm/T 0;
    0 0 0 0 0 1;
    0 Ax*omega_y^2 Ax*2*zeta_y*omega_y 0 -omega_y^2 -2*zeta_y*omega_y;
];

B = [
    0 0;
    0 0;
    Gmx*omega_x^2 0;
    0 0;
    0 0;
    0 Gmy*omega_y^2;
];

C = [
    1 0 0 0 0 0;
    0 0 0 1 0 0;
];

D = [
    0 0;
    0 0
];

n = 6;              % State size
m = 2;              % Input size
p = 2;              % Output size

Ahat = A;
B1hat = [zeros(n,m) B zeros(n,m)];
B2hat = B;
Bhat = [B1hat B2hat];
C1hat = [-C; zeros(m,n); zeros(m,n)];
C2hat = C;
Chat = [C1hat; C2hat];
D11hat = [
    eye(m) -D zeros(m);
    zeros(m) zeros(m) zeros(m);
    eye(m) zeros(m) zeros(m);
];
D12hat = [
     -D;
     eye(m);
     zeros(m);
];
D21hat = [
    zeros(m) D eye(m)
];
D22hat = D;
Dhat = [
     D11hat D12hat;
     D21hat D22hat;
];

m1 = size(B1hat,2);
m2 = size(B2hat,2);
p1 = size(C1hat,1);
p2 = size(C2hat,1);

% sys = ss(Ahat,Bhat,Chat,Dhat);
% sysd = c2d(sys, ts);
% 
% Ad = sysd.A;
% Bd = sysd.B;
% Cd = sysd.C;
% Dd = sysd.D;

save('bench_ABCD','A','B','C','D')