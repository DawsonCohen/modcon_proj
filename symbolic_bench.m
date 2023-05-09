T = 67e-6;

syms Gf Dm zeta_x zeta_y omega_x omega_y Gmx Gmy Ax Ay Vsupply;

A = [
    -1/T Gf*Dm/T 0 0 0 0;
    0 0 1 0 0 0;
    0 -omega_x^2 -2*zeta_x*omega_x 0 Ax*omega_x^2 Ax*2*zeta_x*omega_x;
    0 0 0 -1/T Gf*Dm/T 0;
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

D = zeros(2);

C_ab = [B  A*B  A*A*B];
C_ca = [C; C*A; C*A*A];