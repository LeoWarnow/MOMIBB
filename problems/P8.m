function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,is_convex,is_quadratic] = P8(~)
%P8 A non-convex non-quadratic tri-objective test instance

% Dimension of decision and criterion space
n = 3; % Continuous variables
m = 1; % Integer variables
p = 3; % Dimension criterion space
q = 3; % Number of constraints

% Problem type
is_convex = false;
is_quadratic = false;

% Objective function
f = @(x) [x(1)+x(4);x(2)-x(4);x(3)-exp(x(4))-3];
Df = @(x) [];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = -2.*ones(n+m,1);
ub = 2.*ones(n+m,1);

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) [x(1)^2+x(2)^2-1;x(1)*x(2)*(1-x(3))-1;exp(x(3))-1];
Dg = @(x) [];
end