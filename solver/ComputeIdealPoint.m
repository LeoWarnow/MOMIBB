function [ideal_point,x_res]=ComputeIdealPoint(n,m,p,f,g,Aineq,bineq,Aeq,beq,lb,ub,x0)
%ComputeIdealPoint Computes ideal point

function [c,ceq] = nonlcon_fun(x)
    c = g(x);
    ceq = [];
end
nonlcon = @nonlcon_fun;
options = optimoptions('fmincon','Display','none');

weights = eye(p);
ideal_point = zeros(p,1);
x_res = zeros(n+m,p);
for i=1:p
    w = weights(:,i);
    [x_res(:,i),ideal_point(i)] = fmincon(@(x) w'*f(x),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
end
end