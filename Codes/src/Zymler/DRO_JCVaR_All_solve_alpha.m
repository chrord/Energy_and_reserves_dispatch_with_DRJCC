function[sol] = DRO_JCVaR_All_solve_alpha(si,DRO_param, input, jcc)
%   In this subproblem, the values of x and Y are fixed
%   Solve for alpha

    p = input.p;
    Y = input.Y;
    ru = input.ru;
    rd = input.rd;
    
    x = [p; ru; rd];
    nJCC = size(jcc, 1);
    maxK = 0;
    for j = 1:nJCC
        if size(jcc{j, 1}, 1) > maxK
            maxK = size(jcc{j, 1}, 1);
        end
    end
    
    yalmip('clear')

    % Getting the number of thermals power plants, wind farms, scenarions,
    % transmission lines and nodes
    Nscen = size(si.Wscen,2);

    % Definition of variables
    alpha = sdpvar(nJCC, maxK, 'full');

    s = sdpvar(nJCC, Nscen, 'full');
    lambda = sdpvar(nJCC, 1);     
    tau = sdpvar(nJCC, 1);

    % Constraints set
    CS = [];
    CS = [DRO_param.alpha_min <= alpha <= DRO_param.alpha_max];

    % Run a for-loop to add the constraints related to the joint cvar
    % The set of code below is generic, it can be copied and paste for any
    % structure jcc of interest
    for j = 1:nJCC
        A = jcc{j, 1};
        B = jcc{j, 2};
        C = jcc{j, 3};
        b = jcc{j, 4};
        eps = jcc{j, 5};        
        % find the number of individual chance constraint belonging to this
        K = size(A, 1);        
        CS = [CS, tau(j) <= s(j,:)];
        for k = 1:K
            CS = [CS, (1 - 1/eps)*repmat(tau(j), 1, Nscen) + alpha(j,k)/eps*( repmat(A(k,:)*x - b(k), 1, Nscen) + (B(k,:)*Y+C(k,:))*si.xi ) <= s(j,:) ];
            CS = [CS, norm(alpha(j,k)/eps*(B(k,:)*Y + C(k,:)), DRO_param.dual_norm) <= lambda(j)];
        end
    end
    
    % Optimization options
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0); %, 'verbose', 0

    % Build the objective function 
    Obj = sum(DRO_param.rho*lambda) + (1/Nscen * sum(s(:))); 

    % Solve
    diag = optimize(CS, Obj, optim_options);

    sol.alpha = value(alpha);
    sol.Flag = diag.problem;
    
end