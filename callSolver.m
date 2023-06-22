function [L,U,N,box_struct,it,exitflag,time] = callSolver(problem_name,problem_param,L,U,CUT_MODE,SOL_MODE,EPSILON,OFFSET,plot_result)

%% Load problem
problem = str2func(problem_name);
problem_quadr = [];

%% Handling of incorrect or missing input
[~,~,p,~,f,~,~,~,~,~,~,~,lb,ub,~,~,is_quadratic] = problem(problem_param);

% Applying standard tolerances
if isempty(EPSILON)
    EPSILON = 0.1;
end
if isempty(OFFSET)
    OFFSET = EPSILON*1e-3;
end

% Checking MOMIBB parameters
if isempty(SOL_MODE)
    SOL_MODE = 1;
    disp('Selected SOL_MODE = 1');
elseif SOL_MODE == 2
	if ~is_quadratic
        SOL_MODE = 1;
        warning('No quadratic formulation available. Changing to SOL_MODE = 1.');
    else
        problem_quadr = str2func(problem_name+"_quadr");
        if ~isempty(CUT_MODE) && (CUT_MODE > 0)
            CUT_MODE = 0;
            warning('Cuts are not supported for SOL_MODE = 2. Changing to CUT_MODE = 0.');
        end
	end
end
if isempty(CUT_MODE)
    CUT_MODE = 1;
    disp('Selected cutting strategy using continuous cuts only (CUT_MODE == 1)');
elseif CUT_MODE == 2
    if ~is_quadratic
        CUT_MODE = 1;
        warning('No quadratic formulation available. Changing to CUT_MODE = 1.');
    else
        problem_quadr = str2func(problem_name+"_quadr");
    end
end

% Initialization of L,U
if isempty(L) || isempty(U)
    startbox = infsup(lb, ub);
    B = intval;
    for i=1:p
        B(i) = (1:p==i)*f(startbox);
    end
    if isempty(L)
        L = B.inf';
    end
    if isempty(U)
        U = B.sup';
    end
end
L = L-OFFSET;
U = U+OFFSET;

%% Call MOMIBB
if SOL_MODE == 1
    tic;
    [L,U,N,box_struct,it,exitflag] = MOMIBB(problem,problem_quadr,problem_param,CUT_MODE,L,U,EPSILON,OFFSET);
    time = toc;
elseif SOL_MODE == 2
	tic;
    [L,U,N,box_struct,it,exitflag] = MOMIBB_direct(problem,problem_quadr,problem_param,L,U,EPSILON,OFFSET);
    time = toc; 
end

%% Plot if wanted
if plot_result > 0
    plotBoxes(L,U,p);
    if p < 3
        figure;
        hold on;
        plot(L(1,:),L(2,:),'LineStyle','none','Marker','.','Color','blue');
        plot(U(1,:),U(2,:),'LineStyle','none','Marker','.','Color',[1, 0.4745, 0]);
        grid on;
        xlabel('f_1');
        ylabel('f_2');
    elseif p < 4
        figure;
        hold on;
        plot3(L(1,:),L(2,:),L(3,:),'LineStyle','none','Marker','.','Color','blue');
        plot3(U(1,:),U(2,:),U(3,:),'LineStyle','none','Marker','.','Color',[1, 0.4745, 0]);
        grid on;
        xlabel('f_1');
        ylabel('f_2');
        zlabel('f_3');
    end
end
end