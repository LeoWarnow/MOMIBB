function [sol_x,exitflag]=WSMIQP(n,m,p,q,Aineq,bineq,Aeq,beq,lb,ub,Qfun,qfun,cfun,Qcon,qcon,ccon,x_start,lambda)

clear model;
% Gurobi model
%  minimize
%        lambda(1)(x'Q1x + c1'x) + lambda(2)(x'Q2x + c2'x)
%  subject to
%        g(x) \leq 0
%        x \in [lb,ub]

% Initialize model
model.modelsense = 'min';
model.modelname = 'MixedIntegerHyperplane';

% Variables
model.vtype = [repmat('C',1,n),repmat('I',1,m)];

% Objective function
Qobj = lambda(1)*Qfun{1};
qobj = lambda(1)*qfun{1};
cobj = lambda(1)*cfun{1};
for j=2:p
    Qobj=Qobj+lambda(j)*Qfun{j};
    qobj=qobj+lambda(j)*qfun{j};
    cobj=cobj+lambda(j)*cfun{j};
end
model.Q = sparse(Qobj);
model.obj = qobj;
model.objcon = cobj;

% Constraint functions
% Linear constraints
if isempty(bineq)
    Aineq = zeros(1,m+n);
    bineq = 0;
end
if ~isempty(beq)
    sizeAeq = size(Aeq,1);
    model.A = sparse([Aineq;Aeq]);
    model.rhs = [bineq;beq];
    sizeAineq = size(Aineq,1);
    model.sense = [repmat('<',1,sizeAineq-sizeAeq),repmat('=',1,sizeAeq)];
else
    model.A = sparse(Aineq);
    model.rhs = bineq;
    model.sense = '<';
end
% Quadratic constraints
for j=1:q
    model.quadcon(j).Qc=sparse(Qcon{j});
    model.quadcon(j).q=qcon{j};
    model.quadcon(j).rhs=ccon{j};
    model.quadcon(j).sense='<';
    model.quadcon(j).name=['quadratic_constraint_',num2str(j)];
end

% Upper and lower bounds
model.lb = lb;
model.ub = ub;

% Starting point
model.start = x_start;

% Write model
gurobi_write(model, 'solver/MixedIntegerHyperplane.lp');

% Solve model and return results
clear params;
params.outputflag = 0;
params.NonConvex = 2;
result = gurobi(model, params);
if strcmp(result.status,'OPTIMAL')
    exitflag = 1;
    sol_x = result.x; %TODO: Abgleich mit originalem GurobiCall wegen Dimension
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1;
    sol_x = -1;
else
    exitflag = -2;
    sol_x = -1;
end

end