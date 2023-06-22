function [err]=FeasibilityCheckMIQP(n,m,q,Aineq,bineq,Aeq,beq,lb,ub,Qcon,qcon,ccon,x_start)
%FeasibilityCheckMIQP Checkes if there exists a feasible point within the current box

clear model;

% Initialize model
model.modelsense = 'min';
model.modelname = 'FeasibilityCheckMIQP';

% Variables
model.vtype = [repmat('C',1,n),repmat('I',1,m),'C'];

% Objective function
model.obj = [zeros(1,n+m),1];

% Constraint functions
% Linear constraints
if isempty(bineq)
    Aineq = zeros(1,m+n+1);
    bineq = 0;
else
    size_Aineq = size(Aineq,1);
    Aineq = [Aineq,-ones(size_Aineq,1)];
end
if ~isempty(beq)
    size_Aeq = size(Aeq,1);
    Aeq = [Aeq,zeros(size_Aeq,1)];
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
    model.quadcon(j).Qc=sparse([Qcon{j},zeros(n+m,1);zeros(1,n+m+1)]);
    model.quadcon(j).q=[qcon{j};-1];
    model.quadcon(j).rhs=ccon{j};
    model.quadcon(j).sense='<';
    model.quadcon(j).name=['quadratic_constraint_',num2str(j)];
end

% Upper and lower bounds
model.lb = [lb;-Inf];
model.ub = [ub;Inf];

% Starting point
model.start = [x_start;1];

% Write model
gurobi_write(model, 'solver/FeasibilityCheckMIQP.lp');

% Solve model and return results
clear params;
params.outputflag = 0;
params.NonConvex = 2;
params.timelimit = 3600;
result = gurobi(model, params);
if strcmp(result.status,'OPTIMAL')
    exitflag = 1;
    sol_x = result.x; %TODO: Abgleich mit originalem GurobiCall wegen Dimension
elseif strcmp(result.status,'INFEASIBLE')
    exitflag = -1;
    sol_x = 1;
elseif strcmp(result.status,'TIME_LIMIT')
    exitflag = -3;
    sol_x = 1;
else
    exitflag = -2;
    sol_x = -1;
end
err = sol_x(end);
end