function [ DRO_sol_JCVaR ] = DRO_JCVaR_All( si, DRO_param, jcc )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
    % Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This script implements the sequential programming problem to solve
    % for the joint CVaR
    
    % Notice that this implementation does NOT include the support for xi
    
    % jcc contains all the information for the joint chance constraint of
    % the form P( Ax + (BY+C)xi <= b ) >= 1 - eps
    % For the j-th joint chance constraint:
    % A = jcc{j, 1}
    % B = jcc{j, 2}
    % C = jcc{j, 3}
    % b = jcc{j, 4}
    % eps = jcc{j, 5}
    
    % Initialize the value for alpha
    nJCC = size(jcc, 1);
    maxK = 0;
    for j = 1:nJCC
        if size(jcc{j, 1}, 1) > maxK
            maxK = size(jcc{j, 1}, 1);
        end
    end
    
    % initialize alpha to values of 1
    
    input.alpha = (DRO_param.alpha_min + DRO_param.alpha_max)/2*ones(nJCC, maxK);

    
    iter = 1;

    while iter <= DRO_param.CVaR_max_iter

        % First, fix alpha, solve for x and Y
        sol_xY = DRO_JCVaR_All_solve_xY(si,DRO_param, input, jcc);
        p_sol(:,iter) = sol_xY.p;
        ru_sol(:,iter) = sol_xY.ru;
        rd_sol(:,iter) = sol_xY.rd;
        viol_sol(:,iter) = sol_xY.viol;
        obj(iter) = sol_xY.Obj;

        if iter > 1 && (obj(iter) - obj(iter-1))/obj(iter) > 1e-1
            display('Error: Objective value has to be decreasing\n');
            break
        end

        if iter > 1 && norm((obj(iter) - obj(iter-1))/obj(iter), 2) < DRO_param.tolerance 
            break
        end

        % Update x and Y
        input.p = sol_xY.p;
        input.ru = sol_xY.ru;
        input.rd = sol_xY.rd;
        input.Y = sol_xY.Y;
        
        % solve the alpha subproblem
        sol_alpha = DRO_JCVaR_All_solve_alpha(si,DRO_param, input, jcc);

        % update alpha
        input.alpha = sol_alpha.alpha;
            
        iter = iter + 1;
    end
    
    % Resolve the last time
    DRO_sol_JCVaR = DRO_JCVaR_All_solve_xY(si,DRO_param, input, jcc);
    
    if norm(DRO_sol_JCVaR.viol, 2) < DRO_param.tolerance
        DRO_sol_JCVaR.message = 'Problem is feasible';
        DRO_sol_JCVaR.Flag = 0;
    else
        DRO_sol_JCVaR.message = 'Problem is infeasible';
        DRO_sol_JCVaR.Flag = 12;
    end
end

