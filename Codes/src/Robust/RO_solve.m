function[sol] = RO_solve(si,jcc)

yalmip('clear')

    % Getting the number of thermals power plants, wind farms, scenarions,
    % transmission lines and nodes
    Nunits = size(si.Pmax,1);
    Nwind = size(si.Wmax,1);
    Nscen = size(si.Wscen_RO,2);
   

    % Definition of variables
    p = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    ru = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    rd = sdpvar(Nunits, 1); % Day-ahead power production from thermal power plants
    Y = sdpvar(Nunits, Nwind, 'full'); % Linear decision rule for real-time power production
    beta = sdpvar(1,1); % Worst case cost
    
    % create x by stacking up p, ru and rd
    x = [p; ru; rd];
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    CS = [CS, si.Pmin <= p - rd, p + ru <= si.Pmax, 0 <= ru <= si.ResCap, 0 <= rd <= si.ResCap];
    CS = [CS, sum(p) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(Y, 1) == -si.Wmax'];
    
    % Run a for-loop to add the constraints related to the individual cvar
    % The set of code below is generic, it can be copied and paste for any
    % structure jcc of interest
    
    % find the number of Individual chance constraints we have
    nICC = 0;
    for j=1:size(jcc, 1)
        nICC = nICC + size(jcc{j, 1}, 1);
    end
    
    for j=1:size(jcc, 1)
        A_C{j,1} = jcc{j,1};
        B_C{j,1} = jcc{j,2};
        C_C{j,1} = jcc{j,3};
        b_C{j,1} = jcc{j,4};
    end
    A = cell2mat(A_C);
    B = cell2mat(B_C);
    C = cell2mat(C_C);
    b = cell2mat(b_C);

    
    for j = 1:nICC
    CS = [CS, repmat(A(j,:)*x - b(j), 1, Nscen) + (B(j,:)*Y+C(j,:))*si.zi  <= 0 ];
    end
    
    % Build the objective function 
    Obj = si.Cr1'*ru + si.Cr2'*rd + si.C'*p + beta; 
    CS = [CS, si.C'*Y*si.zi <= beta];

    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.p = value(p);
    sol.Y = value(Y);
    sol.ru = value(ru);
    sol.rd = value(rd);
    sol.Y = value(Y);
    sol.fy = si.Qg*value(p) + si.Qw*si.DiagWmax * si.mu - si.Qd*si.D;
    sol.fY = si.Qg*value(Y) + si.Qw*si.DiagWmax;
    sol.q = si.PG * value(p);
    sol.qY = si.PG * value(Y);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end