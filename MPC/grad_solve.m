%% Gradient Descent
function [x] = grad_solve(F, gradF, x0, alpha, beta)
    arguments
        F
        gradF
        x0 (:,1)
        alpha double {mustBeNonempty} = 0.1
        beta double {mustBeNonempty} = 0.7
    end

    eta = 1e-7;
    x = x0;
    while(true)
        Dx = -gradF(x);
        if(norm(Dx) < eta); break; end % Stopping criterion
        t = BTLS(F,x,Dx,alpha, beta);
        x = x + t * Dx;
    end
end

%% Line Search
function t = BTLS(F,x,Dx,alpha,beta)
    t = 1;
    f = F(x);
    while(F(x+t*Dx) > f-alpha*t*(Dx'*Dx))
        t = beta*t;
    end
end