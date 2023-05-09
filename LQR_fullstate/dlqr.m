function [K,P] = dlqr(A,B,Q,R)

P = idare(A,B,Q,R,[],[]);
K = -(R+B'*P*B)\B'*P*A;

end