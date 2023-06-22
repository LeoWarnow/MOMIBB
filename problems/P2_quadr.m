function [Qfun,qfun,cfun,Qcon,qcon,ccon] = P2_quadr(params)

% Dimension of decision and criterion space
n = params(1); % Continuous variables
m = params(2); % Integer variables
r = n+m;
assert(mod(n,2)==0,'Number of continuous variables has to be even.')
assert(mod(n,2)==0,'Number of continuous variables has to be even.')

% Objective function
Qfun{1} = zeros(r);
qfun{1} = [ones(1,n/2),zeros(1,n/2),ones(1,m/2),zeros(1,m/2)]';
cfun{1} = 0;
Qfun{2} = zeros(r);
qfun{2} = [zeros(1,n/2),ones(1,n/2),zeros(1,m/2),ones(1,m/2)]';
cfun{2} = 0;

% Constraints
Qcon{1} = diag([-ones(1,n),zeros(1,m)]);
qcon{1} = zeros(r,1);
ccon{1} = -1;
Qcon{2} = diag([zeros(1,n),ones(1,m)]);
qcon{2} = zeros(r,1);
ccon{2} = 9;
end
