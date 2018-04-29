%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script for the robust optimization model
%

clear all; close all; clc;

% call startup to add the necessary path
startup;

tic

%%

% Data input

RTS_Data2;

% Number of individual runs (number of coupled datasets in the numerical
% study)

IR_max = 100;
IR_sim = 100;

% Number of out of sample data for each individual run (N') for testing
% dataset

OOS_max = 200;
OOS_sim = 100;

% Number of maximum sample size (N)

N_max = 1000;

% Number of sample data in training dataset (N)

N = 100;

% Total number of data 

Nscen = IR_max * (N_max + OOS_max);

% Generation of data 

rng(4,'twister')

% Number of wind farms

wf=[1:6];

% Loading the historical data for wind farms

wff=AV_AEMO2(:,wf);

% Cutting off very extreme values

cut_off_eps = 1e-2;
wff(wff<cut_off_eps) = cut_off_eps;
wff(wff>(1-cut_off_eps)) = 1 - cut_off_eps;

% Logit-normal transformation (Eq. (1) in ref. [31])

yy=log(wff./(1-wff));

% Calculation of mean and variance, note that we increase the mean to have
% higher wind penetration in our test-case

mu = mean(yy)+1.5;
sigma_m=cov(yy);
sigma_m=sigma_m./(std(yy)'*std(yy));

% Inverse of logit-normal transformation (Eq. (2) in ref. [31])

R = chol(sigma_m);
y = repmat(mu,Nscen,1) + randn(Nscen,size(WindDATA,1))*R;
Wind = (1+exp(-y)).^(-1);

% Checking correlation, mean and true mean of data

corrcoef(Wind);
mean(Wind);
true_mean_Wind = (1+exp(-mu)).^(-1);

% Reshaping the data structure

nWind = Wind';
nWind = reshape(nWind,size(WindDATA,1), N_max+OOS_max, IR_max);

% Initializing the matrices to gather final results

RO_TC = NaN(IR_sim);

for j = 1:IR_sim
    display('out of sample iteration:');
    j
  
    % For each coupled dataset, we pick N and N' samples
    WPf_max = nWind(:,1:N_max,j)';
    WPr_max = nWind(:,N_max+1:N_max+OOS_max,j)';
    WPf = WPf_max(1:N,:);
    WPr = WPr_max(1:OOS_sim,:);
    
    % Build the corresponding data related to wind power production
    all = [1:N];
    system_info.Wscen = WPf(all,:)';
    system_info.mu = mean(system_info.Wscen,2); 
    system_info.xi = system_info.Wscen - repmat(system_info.mu, 1, size(system_info.Wscen,2));

    % Build the corresponding data for RO
    system_info.Wscen_RO = de2bi(0:2^size(system_info.Wmax,1)-1)';
    system_info.Wexp_RO = mean(system_info.Wscen_RO,2); 
    system_info.zi = system_info.Wscen_RO - repmat(system_info.mu, 1, size(system_info.Wscen_RO,2));
    
    % Calculation of A,B,C,b matrices for joint chance constraints
    CC_jcc = CC_matrices(system_info, DRO_param);
    jcc = CC_jcc.jcc;
    
    tic
    
    % Solve robust optimization model
    
    RO_sol = RO_solve(system_info, jcc);
    RO_p_DA{j} = RO_sol.p;
    RO_ru{j} = RO_sol.ru;
    RO_rd{j} = RO_sol.rd;
    RO_obj{j} = RO_sol.Obj;
    RO_flag{j} = RO_sol.Flag;
    RO_Y{j} = RO_sol.Y * system_info.zi;
    RO_Qy{j} = RO_sol.q;
    RO_QY{j} = RO_sol.qY * system_info.zi;
    RO_Fy{j} = RO_sol.fy;
    RO_FY{j} = RO_sol.fY * system_info.zi;
    
    TimeRO(j) = toc;
    
        tic
        % Loop for each out-of-sample realization 
        for k = 1:OOS_sim
            system_info.Wreal = WPr(k,:)';
            system_info.DWreal = system_info.Wreal - system_info.mu;

            % Solve real-time optimal power flow
            RT_solution_CVaR = RT_solve_R(system_info,RO_sol.p,RO_sol.ru,RO_sol.rd);
            RT_Obj_IR(j,k) = RT_solution_CVaR.Obj_RT;  
            RT_lshed{j,k} = RT_solution_CVaR.lshed_RT;
            RT_flow{j,k} = RO_sol.fy + RO_sol.fY * system_info.DWreal;
            RT_p{j,k} = RO_sol.Y * system_info.DWreal;
            RT_q{j,k} = RO_sol.q + RO_sol.qY * system_info.DWreal;
            RT_flag(j,k) = RT_solution_CVaR.Flag;
           
        end
        
        TimeOOS(j) = toc;
        
        % Calculation of expected cost
        
        if RO_sol.Flag == 0
            RO_TC(j) = system_info.Cr1'*RO_sol.ru + system_info.Cr2'*RO_sol.rd + mean(RT_Obj_IR(j,:));
        else
            RO_TC(j) = NaN;
        end
       

 
end


Time = toc
