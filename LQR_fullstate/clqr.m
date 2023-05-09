function [K,P] = clqr(A,B,Q,R)

P = care(A,B,Q,R,[],[]);
K = -R\(B'*P);
% K = -(R+B'*P*B)\B'*P*A;

end