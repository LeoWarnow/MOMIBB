function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,is_convex,is_quadratic] = P3(params)
%P3 A simple non convex example

% Dimension of decision and criterion space
n = params(1); % Continuous variables
m = 1; % Integer variables
p = 2; % Dimension criterion space
q = 1; % Number of constraints
assert(mod(n,2)==0,'Number of continuous variables has to be even.')

% Problem type
is_convex = false;
is_quadratic = false;

% Objective function
f = @(x) [sum(x(1:n/2))+x(n+m);sum(x(n/2+1:n))-exp(x(n+m))];
Df = @(x) [ones(1,n/2),zeros(1,n/2),1;zeros(1,n/2),ones(1,n/2),-exp(x(n+m))];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = [zeros(n,1);-4.*ones(m,1)];
ub = 1.*ones(n+m,1);

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) [-sum(x(1:n).^2)+1];
Dg = @(x) [-2.*x(1:n)',0];
end