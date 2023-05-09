%% Newton's Method
function [x, hist] = newton_solve(F,gradF,hessF,val,x0,eps,alpha,beta)
    arguments
        F
        gradF
        hessF
        val
        x0 (:,1)
        eps double = 1e-4
        alpha double = 0.01
        beta double = 0.5
    end
    
    x = x0;
    hist = [];
    while(true)
        % Update Variables
        Grad = gradF(x);
        Hess = hessF(x);
        % Save History
        hist = [hist; x'];
        [Dx, lam_squared] = newton_step(Grad, Hess);
        if(lam_squared/2 < eps); break; end
        t = BTLS(F,gradF,val,x,Dx,alpha,beta);
        x = x + t*Dx;
    end
end

%% Newton Step
function [Dx, lam_squared] = newton_step(Grad, Hess)
    Dx = -Hess\Grad;
    lam_squared = -Grad' * Dx;
end

%% Line Search
function t = BTLS(F,gradF,val,x,Dx,alpha,beta)
    t = 1;
    f = F(x);
    lambda = -gradF(x)'*Dx;
    
    while (~val(x+t*Dx))
        t = beta*t;
    end

    while(F(x+t*Dx) > f+alpha*t*lambda)
        t = beta*t;
    end
end