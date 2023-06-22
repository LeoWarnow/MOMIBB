function [ideal_point,x_res]=ComputeIdealPointMIQP(n,m,p,q,Aineq,bineq,Aeq,beq,lb,ub,Qfun,qfun,cfun,Qcon,qcon,ccon,x_start)
%ComputeIdealPointMIQP Computes ideal point

clear model;

% Initialize model
model.modelsense = 'min';
model.modelname = 'MixedIntegerHyperplane';

% Variables
model.vtype = [repmat('C',1,n),repmat('I',1,m)];

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

% Solver loop
weights = eye(p);
ideal_point = zeros(p,1);
exitflag = zeros(p,1);
x_res = zeros(n+m,p);
for i=1:p    
    % Objective function
    w = weights(:,i);
    Qobj = w(1)*Qfun{1};
    qobj = w(1)*qfun{1};
    cobj = w(1)*cfun{1};
    for j=2:p
        Qobj=Qobj+w(j)*Qfun{j};
        qobj=qobj+w(j)*qfun{j};
        cobj=cobj+w(j)*cfun{j};
    end
    model.Q = sparse(Qobj);
    model.obj = qobj;
    model.objcon = cobj;
    
    % Write model
    gurobi_write(model, 'solver/IdealPointMIQP.lp');

    % Solve model and return results
    clear params;
    params.outputflag = 0;
    params.NonConvex = 2;
    result = gurobi(model, params);
    if strcmp(result.status,'OPTIMAL')
        exitflag(i) = 1;
        x_res(:,i) = result.x; %TODO: Abgleich mit originalem GurobiCall wegen Dimension
        ideal_point(i) = result.objval;
    elseif strcmp(result.status,'INFEASIBLE')
        exitflag(i) = -1;
        x_res(:,i) = -1;
    else
        exitflag(i) = -2;
        x_res(:,i) = -1;
    end
end
end