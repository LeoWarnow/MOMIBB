function [L,U,N,box_struct,it,flag_main] = MOMIBB(problem,problem_quadr,problem_param,CUT_MODE,L,U,EPSILON,OFFSET)

%% Initialization phase
[n,m,p,q,f,g,~,~,Aineq,bineq,Aeq,beq,lb,ub,~,is_convex,~] = problem(problem_param);
f_orig = f;
g_orig = g;
if isempty(Aineq)
    is_feasible = @(x) all(g_orig(x) <= 0);
else
    is_feasible = @(x) all([g_orig(x) <= 0;Aineq*x <= bineq]);
end
box_struct = cell(0,4); %(X,LB,ideal,finished)
N = [];
it = 0;
TIME_LIMIT = 3600;
flag_main = -2;

% Load quadratic formulation if needed
if CUT_MODE == 2
    [Qfun,qfun,cfun,Qcon,qcon,ccon] = problem_quadr(problem_param);
end

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
        
        % For nonconvex problems use aBB
        if ~is_convex
            x_aBB = hessianinit(current_subbox);
            psi = @(x) (1/2)*((subbox_lb - x)'*(subbox_ub - x));
            alpha_f = zeros(p,1); 
            alpha_g = zeros(q,1);
            for j=1:p
                hess = (1:p==j)*f_orig(x_aBB);
                A = hess.hx;
                A_lb = A.inf;
                A_ub = A.sup;
                temp = max(abs(A_lb),abs(A_ub));
                temp(1:(n+m)+1:(n+m)*(n+m)) = 0;
                beta_j = diag(A_lb) - sum(temp,2);
                beta = min(beta_j);    
                alpha_f(j) = max([0 (-1)*beta]);
            end
            for j=1:q
                hess = (1:q==j)*g_orig(x_aBB);
                if hess == 0
                    alpha_g(j) = 0;
                    continue;
                end
                A = hess.hx;
                A_lb = A.inf;
                A_ub = A.sup;
                temp = max(abs(A_lb),abs(A_ub));
                temp(1:(n+m)+1:(n+m)*(n+m)) = 0;
                beta_j = diag(A_lb) - sum(temp,2);
                beta = min(beta_j);    
                alpha_g(j) = max([0 (-1)*beta]);
            end
            f = @(x) f_orig(x) + alpha_f*psi(x);
            g = @(x) g_orig(x) + alpha_g*psi(x);
        end

        if all(subbox_lb==subbox_ub)
            % Handling of boxes containing a single point
            subbox_x = subbox_lb;
            subbox_ideal = subbox_y;
            subbox_L = subbox_ideal;
            if is_feasible(subbox_x)
                subbox_y = f_orig(subbox_x);
                [~,flag] = updateNDS(N,subbox_y);
                if flag
                    N = [N subbox_y];
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
            err = FeasibilityCheck(q,g,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,subbox_x0);
            if err > OFFSET
                subbox_discard = true;
            else
                [subbox_ideal,eff_candidates] = ComputeIdealPoint(n,m,p,f,g,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,subbox_x0);
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
                const_cut = diag(H(1:p,:)'*H(p+1:end-1,:));
                depth_cut = min(H(1:p,:)'*U - const_cut);
                ub_index = (depth_cut>-OFFSET);
                U_loop = U(:,ub_index);
                depth_cut = depth_cut(ub_index);
                size_U_loop = size(U_loop,2);
                for l=1:size_U_loop
                    [depth,u_index] = max(depth_cut);
                    if (toc > TIME_LIMIT), break; end %recognize time limit early
                    if depth >= 0
                        u = U_loop(:,u_index);
                        [x_cut,t_cut,~,lagrange_multi] = PS(f,g,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,subbox_x0,u,ones(p,1));
                        lambda = lagrange_multi.ineqnonlin(1:p);
                        lambda = lambda/norm(lambda);
                        x_cut(n+1:n+m) = round(x_cut(n+1:n+m));
                        z_cut = u + t_cut.*ones(p,1);
                        if is_feasible(x_cut)
                            y_cut = f_orig(x_cut);
                            N = [N y_cut];
                            U = updateLUB3(U,y_cut);
                        end
                        if CUT_MODE == 0 || any(isnan(lambda))
                            depth_cut(u_index) = -1;
                            if t_cut < OFFSET
                                subbox_discard = false;
                                break;
                            end
                        elseif t_cut < OFFSET && (CUT_MODE == 2)
                            [x_cut,exitflag] = WSMIQP(n,m,p,q,Aineq,bineq,Aeq,beq,subbox_lb,subbox_ub,Qfun,qfun,cfun,Qcon,qcon,ccon,subbox_x0,lambda);
                            if exitflag == -1
                                subbox_discard = true;
                                break;
                            elseif exitflag == 1
                                y_cut = f_orig(x_cut);
                                z_cut = y_cut;
                                if is_feasible(x_cut)
                                    N = [N y_cut];
                                    U = updateLUB3(U,y_cut);
                                end
                                H = [H,[lambda;z_cut;3]];
                                const_cut(end+1) = lambda'*z_cut;
                                depth_cut = min([depth_cut;lambda'*U_loop-lambda'*z_cut]);
                                if lambda'*u >= lambda'*z_cut
                                    subbox_discard = false;
                                    break;
                                end
                            else
                                subbox_discard = false;
                                break;
                            end                                
                        elseif t_cut < OFFSET
                            H = [H,[lambda;z_cut;2]];
                            subbox_discard = false;
                            break;
                        else
                            H = [H,[lambda;z_cut;2]];
                            const_cut(end+1) = lambda'*z_cut;
                            depth_cut = min([depth_cut;lambda'*U_loop-lambda'*z_cut]);
                        end
                    else
                        subbox_discard = true;
                        break;
                    end
                end
            end      
        end

        % Update the box structure
        if ~subbox_discard
            box_struct(end+1,:) = {current_subbox,H,subbox_L,true}; %(X,LB,ideal,finished)
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