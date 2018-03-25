function [ RT_sol ] = RT_solve_R(si, x, ru, rd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given x (the day-ahead commitment), solve the real time dispatch

yalmip('clear')

% Getting the number of thermals power plants, wind farms, scenarions,
% transmission lines and nodes
Nunits = size(si.Pmax,1);
Nwind = size(si.Wmax,1);
Nlines = size(si.F,1);
Nnodes = size(si.AG,1);
Nloads = size(si.D,1);

% Definition of variables
dp = sdpvar(Nunits, 1); % Real-time power production from thermal power plants
lshed = sdpvar(Nloads, 1); % Real-time load shedding
wsp = sdpvar(Nwind, 1);  % Real-time wind spilling
fi_real = sdpvar(Nnodes, 1); % Day-ahead injection at each node

% Constraints set
CS = [];
CS = [CS, si.Pmin <= x + dp <= si.Pmax, -rd <= dp <= ru, -si.F <= si.PTDF*fi_real <= si.F, sum(fi_real) == 0];
CS = [CS, si.AG*(x + dp) + si.AW * si.DiagWmax*si.Wreal - si.AD*(si.D - lshed) - si.AW * wsp == fi_real, 0 <= lshed <= si.D];
CS = [CS, 0 <= si.AW * wsp <= si.AW * si.DiagWmax*si.Wreal];
CS = [CS, 0 <= si.PG * (x+dp) <= si.FP];


% Build the onjective function 
Obj_real = si.C' * (x+dp) + si.Clsh'*lshed; 

% Optimization options
optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',500,'verbose',0);

% Solve
sol = optimize(CS, Obj_real, optim_options);

RT_sol.p_RT = x + value(dp);
RT_sol.lshed_RT = value(lshed);
RT_sol.wsp_RT = value(wsp);
RT_sol.flow_RT = si.PTDF*value(fi_real);
RT_sol.Obj_RT = value(Obj_real);
RT_sol.Flag = sol.problem;

end

