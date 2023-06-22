function [err]=FeasibilityCheck(q,g,Aineq,bineq,Aeq,beq,lb,ub,x_start)
%FeasibilityCheck Checkes if there exists a feasible point within the current box

% If there are any constraints for the original problem then they nned to
% be extended inluding the new variable
if ~isempty(Aineq)
    size_Aineq = size(Aineq,1);
    Aineq = [Aineq,-ones(size_Aineq,1)];
end
if ~isempty(Aeq)
    size_Aeq = size(Aeq,1);
    Aeq = [Aeq,zeros(size_Aeq,1)];
end

lb = [lb;-Inf];
ub = [ub;Inf];

function [c,ceq] = nonlcon_fun(x)
    c = g(x(1:end-1)) - x(end).*ones(q,1);
    ceq = [];
end
nonlcon = @nonlcon_fun;
options = optimoptions('fmincon','Display','none');

x0 = [x_start;1];

solution = fmincon(@(x) x(end),x0,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,options);
err = solution(end);
end