function [L,U,N,box_struct,it,flag_main] = MOMIBB_direct(problem,problem_quadr,problem_param,L,U,EPSILON,OFFSET)

%% Initialization phase
[n,m,p,q,f,g,~,~,Aineq,bineq,Aeq,beq,lb,ub,~,~,~] = problem(problem_param);
[Qfun,qfun,cfun,Qcon,qcon,ccon] = problem_quadr(problem_param);
f_orig = f;
g_orig = g;
if isempty(Aineq)
    is_feasible = @(x) all(g_orig(x) <= 0);
else
    is_feasible = @(x) all([g_orig(x) <= 0;Aineq*x <= bineq]);
end
box_struct = cell(0,4); %(X,LB,ideal,x_feas)
N = [];
it = 0;
TIME_LIMIT = 3600;
flag_main = -2;

%% Computation for the first box
current_box = infsup(lb, ub);
box_ideal = min(L,[],2);
box_struct(1,:) = {current_box,[eye(p);repmat(box_ideal,1,p);ones(1,p)],L,true};
size_box_struct = 1;

%% Main loop
tic;
while (toc<TIME_LIMIT)
    % Increase number of iterations
    it = it+1;

    % Select box with largest shortest edge
    s = zeros(size_box_struct,1);
    for i = 1:size_box_struct
        s_temp = 0;
        for l = box_struct{i,3}
            U_temp_index = all(l<(U-EPSILON+OFFSET));
            if any(U_temp_index)
                s_temp = max(min(U(:,U_temp_index)-l));
            end
            s(i) = max([s_temp,s(i)]);
        end
    end
    [enclosure_width,current_box_index] = max(s);
    current_box = box_struct{current_box_index,1};

    if enclosure_width < EPSILON
        flag_main = 1;
        break;
    end

    % Branching step
    current_box_diam = diam(current_box);
    [~,branching_index] = max(current_box_diam(end:-1:1));
    branching_index = n+m+1-branching_index;
    subboxes = {current_box,current_box};
    if branching_index < n+1
        subboxes{1}(branching_index) = infsup(inf(current_box(branching_index)),mid(current_box(branching_index)));
        subboxes{2}(branching_index) = infsup(mid(current_box(branching_index)),sup(current_box(branching_index)));
    elseif mod(inf(current_box(branching_index))+sup(current_box(branching_index)),2)
        subboxes{1}(branching_index) = infsup(ceil(inf(current_box(branching_index))),floor(mid(current_box(branching_index))));
        subboxes{2}(branching_index) = infsup(ceil(mid(current_box(branching_index))),floor(sup(current_box(branching_index))));
    else
        subboxes{1}(branching_index) = infsup(ceil(inf(current_box(branching_index))),floor(mid(current_box(branching_index))));
        subboxes{2}(branching_index) = infsup(ceil(mid(current_box(branching_index)))+1,floor(sup(current_box(branching_index))));
    end
    
    % Remove old box
    box_struct(current_box_index,:) = [];
    
    % Bounds & Cuts
    for i=1:2
        current_subbox = subboxes{i};
        subbox_x0 = mid(current_box);  
        subbox_lb = inf(current_box);
        subbox_ub = sup(current_box);

        if all(subbox_lb==subbox_ub)
            % Handling of boxes containing a single point
            subbox_x = subbox_lb;
            subbox_ideal = subbox_y;
            subbox_L = subbox_ideal;
            if is_feasible(subbox_x)
                subbox_y = f_orig(subbox_x);
                [indexlist,flag] = updateNDS(N,subbox_y);
                if flag
                    N = [N subbox_y];
                    N = N(:,indexlist);
                    U = updateLUB3(U,subbox_y);
                    subbox_discard = false;
                else
                    subbox_discard = true;
                end
            else
                subbox_discard = true;
            end
        else
            % Handling of non-singleton boxes
            subbox_ideal = zeros(p,1);
            err = FeasibilityCheckMIQP(n,m,q,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,Qcon,qcon,ccon,subbox_x0);
            if err > OFFSET
                subbox_discard = true;
            else
                [subbox_ideal,eff_candidates] = ComputeIdealPointMIQP(n,m,p,q,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,Qfun,qfun,cfun,Qcon,qcon,ccon,subbox_x0);
                eff_candidates = [eff_candidates(1:n,:);round(eff_candidates(n+1:end,:))];
                for subbox_x = eff_candidates
                    if is_feasible(subbox_x)
                        subbox_y = f_orig(subbox_x);
                        N = [N subbox_y];
                        U = updateLUB3(U,subbox_y);
                    end
                end
                H=[eye(p);repmat(subbox_ideal,1,p);ones(1,p)];
                subbox_discard = true;
                subbox_L = subbox_ideal;
                [dist,idx] = sort(min(U - subbox_ideal),'descend');
                ub_index = idx(dist>-OFFSET);
                for u=U(:,ub_index)
                    if (toc > TIME_LIMIT), subbox_discard = false; break; end %recognize time limit early
                    [x_cut,t_cut,exitflag] = PSMIQP(n,m,p,q,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,Qfun,qfun,cfun,Qcon,qcon,ccon,subbox_x0,u,ones(p,1));
                    x_cut(n+1:n+m) = round(x_cut(n+1:n+m));
                    if exitflag > -0.5
                        y_cut = f_orig(x_cut);
                        N = [N y_cut];
                        U = updateLUB3(U,y_cut);
                    end
                    if t_cut < OFFSET
                        subbox_discard = false;
                        break;
                    end
                end
            end      
        end

        % Update the box structure
        if ~subbox_discard
            box_struct(end+1,:) = {current_subbox,H,subbox_L,true}; %(X,LB,ideal,x_feas)
        end
    end 

    size_box_struct = size(box_struct,1);
    if size_box_struct < 1
        flag_main = 0;
        break;
    end
end
if toc > TIME_LIMIT
    flag_main = -1;
end

%% Compute lower bound set
L = [box_struct{:,3}];
L_indexlist = computeNDS(L);
L = L(:,L_indexlist);
end