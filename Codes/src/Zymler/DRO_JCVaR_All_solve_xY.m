function[sol] = DRO_JCVaR_All_solve_xY(si, DRO_param, input, jcc)
%   In this subproblem, the values of alpha are fixed
%   Solve for x and Y
%   x contains 3 variables: p, r_u, r_d

    alpha = input.alpha;
    
    
    yalmip('clear')

    % Getting the number of thermals power plants, wind farms, scenarions,
    % transmission lines and nodes
    Nunits = size(si.Pmax,1);
    Nwind = size(si.Wmax,1);
    Nscen = size(si.Wscen,2);

    % Definition of variables
    p = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    ru = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    rd = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    %y = sdpvar(Nunits, 1); % Linear decision rule for real-time power production 
    Y = sdpvar(Nunits, Nwind, 'full'); % Linear decision rule for real-time power production
    
    
    
    s_obj = sdpvar(1, Nscen); % sigma variable for obj
    lambda_obj = sdpvar(1, 1); % lambda variable for obj

    
    
    
    % create x by stacking up p, ru and rd
    x = [p; ru; rd];
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    CS = [CS, si.Pmin <= p - rd, p + ru <= si.Pmax, 0 <= ru <= si.ResCap, 0 <= rd <= si.ResCap];
    CS = [CS, sum(p) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(Y, 1) == -si.Wmax'];
    
    % Run a for-loop to add the constraints related to the joint cvar
    % The set of code below is generic, it can be copied and paste for any
    % structure jcc of interest
    
    % find the number of Joint chance constraints we have
    nJCC = size(jcc, 1);
    % create variables
    s = sdpvar(nJCC, Nscen, 'full');
    lambda = sdpvar(nJCC, 1);     
    tau = sdpvar(nJCC, 1); 
    viol = sdpvar(nJCC, 1);
    
    CS = [CS, viol >= 0];
    for j = 1:nJCC
        A = jcc{j, 1};
        B = jcc{j, 2};
        C = jcc{j, 3};
        b = jcc{j, 4};
        eps = jcc{j, 5};
        
        % find the number of individual chance constraint belonging to this
        K = size(A, 1);
        
        CS = [CS, DRO_param.rho*lambda(j) + sum(s(j, :))/Nscen <= viol(j)];
        CS = [CS, tau(j) <= s(j,:)];
        for k = 1:K
            CS = [CS, (1 - 1/eps)*repmat(tau(j), 1, Nscen) + alpha(j,k)/eps*( repmat(A(k,:)*x - b(k), 1, Nscen) + (B(k,:)*Y+C(k,:))*si.xi ) <= s(j,:) ];
            CS = [CS, norm(alpha(j,k)/eps*(B(k,:)*Y + C(k,:)), DRO_param.dual_norm) <= lambda(j)];
        end
    end
   
    % Build the objective function 
    Obj = si.Cr1'*ru + si.Cr2'*rd + si.C'*p + DRO_param.rho*lambda_obj + 1/Nscen * sum(s_obj) + DRO_param.penalty*sum(si.C)*sum(viol); 
    CS = [CS, si.C'*Y*si.xi <= s_obj];
    CS = [CS, norm( Y' * si.C, DRO_param.dual_norm) <= lambda_obj];
    
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.p = value(p);
    sol.Y = value(Y);
    sol.ru = value(ru);
    sol.rd = value(rd);
    sol.viol = value(viol);
    sol.Y = value(Y);
    sol.fy = si.Qg*value(p) + si.Qw*si.DiagWmax * si.mu - si.Qd*si.D;
    sol.fY = si.Qg*value(Y) + si.Qw*si.DiagWmax;
    sol.q = si.PG * value(p);
    sol.qY = si.PG * value(Y);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;

end