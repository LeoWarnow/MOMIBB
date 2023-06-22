function [n,m,p,q,f,g,Df,Dg,Aineq,bineq,Aeq,beq,lb,ub,x0,is_convex,is_quadratic] = P2(params)
%P2 A simple non convex example

% Dimension of decision and criterion space
n = params(1); % Continuous variables
m = params(2); % Integer variables
p = 2; % Dimension criterion space
q = 2; % Number of constraints
assert(mod(n,2)==0,'Number of continuous variables has to be even.')
assert(mod(m,2)==0,'Number of integer variables has to be even.')

% Problem type
is_convex = false;
is_quadratic = true;

% Objective function
f = @(x) [sum(x(1:n/2))+sum(x(n+1:n+m/2));sum(x(n/2+1:n))+sum(x(n+m/2+1:n+m))];
Df = @(x) [ones(1,n/2),zeros(1,n/2),ones(1,m/2),zeros(1,m/2);zeros(1,n/2),ones(1,n/2),zeros(1,m/2),ones(1,m/2)];

% Linear constraints (Aineq*x <= bineq, Aeq*x = beq)
Aineq = [];
bineq = [];
Aeq = [];
beq = [];

% Lower and upper bounds (lb <= x <= ub)
lb = [zeros(n,1);-3.*ones(m,1)];
ub = 3.*ones(n+m,1);

% Start point x0
x0 = ceil((lb+ub)/2);

% Non-linear constraints (g(x) <= 0)
g = @(x) [-sum(x(1:n).^2)+1;sum(x(n+1:n+m).^2)-9];
Dg = @(x) [-2.*x(1:n)',zeros(1,m);zeros(1,n),2.*x(n+1:n+m)'];
end