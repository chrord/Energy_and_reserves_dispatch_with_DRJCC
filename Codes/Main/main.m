%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script
%

clear all; close all; clc;

% call startup to add the necessary path
startup;

tic

%%

% Data input

RTS_Data2;

% DRO input

% Definition of dual norm used in the Wasserstein metric, value of \epsilon
% and the parameters for the algorithm in the Zymler approximation

DRO_param.dual_norm = 'inf'; % dual norm
DRO_param.eps_joint_cvar = 0.05; % \epsilon
DRO_param.CVaR_max_iter = 40; % MaxIter
DRO_param.tolerance = 1e-1; % \eta
DRO_param.alpha_max = 1000; % \overline{\delta}
DRO_param.alpha_min = 1e-4; % \undeline{\delta}
DRO_param.penalty = 1e6; % BigM

% Values of rho, the finite positive-valued vector P is defined

% Vector for Bonferroni approximation

rho_vectorC = [0 linspace(0.0001, 0.0025, 24)];

% Vector for Zymler approximation

rho_vectorJC = [0 0.0001 linspace(0.001, 0.01, 23)]; 

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

Joint_CVaR_Obj_IR = zeros(IR_sim, OOS_sim, length(rho_vectorC));
CVaR_Obj_IR = zeros(IR_sim, OOS_sim, length(rho_vectorJC));
ICC_TC = NaN(IR_sim,length(rho_vectorC));
JCC_TC = NaN(IR_sim,length(rho_vectorJC));

% Loop for each individual run for 100 coupled datasets

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
    
    % Calculation of A,B,C,b matrices for joint chance constraints
    CC_jcc = CC_matrices(system_info, DRO_param);
    jcc = CC_jcc.jcc;
    
    % Loop for each value of \rho in P vector
    for i = 1:length(rho_vectorC) 
             
        % optimize for each value of rho for Bonferroni approximation
        
        DRO_param.rho = rho_vectorC(i);
          
        DRO_ICC_CVaR = DRO_CVaR_ICC(system_info, DRO_param, jcc);
        ICC_p_DA{j, i} = DRO_ICC_CVaR.p;
        ICC_ru{j, i} = DRO_ICC_CVaR.ru;
        ICC_rd{j, i} = DRO_ICC_CVaR.rd;
        ICC_obj{j, i} = DRO_ICC_CVaR.Obj;
        ICC_flag{j, i} = DRO_ICC_CVaR.Flag;
        CVaR_Y{j,i} = DRO_ICC_CVaR.Y * system_info.xi;
        CVaR_Qy{j,i} = DRO_ICC_CVaR.q;
        CVaR_QY{j,i} = DRO_ICC_CVaR.qY * system_info.xi;
        CVaR_Fy{j,i} = DRO_ICC_CVaR.fy;
        CVaR_FY{j,i} = DRO_ICC_CVaR.fY * system_info.xi;
        
        % optimize for each value of rho for Zymler approximation
        
        DRO_param.rho = rho_vectorJC(i);
        
        DRO_JCC_CVaR = DRO_JCVaR_All(system_info, DRO_param, jcc);
        JCC_p_DA{j, i} = DRO_JCC_CVaR.p;
        JCC_ru{j, i} = DRO_JCC_CVaR.ru;
        JCC_rd{j, i} = DRO_JCC_CVaR.rd;
        JCC_obj{j, i} = DRO_JCC_CVaR.Obj;
        JCC_flag{j, i} = DRO_JCC_CVaR.Flag;            
        Joint_CVaR_Y{j,i} = DRO_JCC_CVaR.Y * system_info.xi;
        Joint_CVaR_Qy{j,i} = DRO_JCC_CVaR.q;
        Joint_CVaR_QY{j,i} = DRO_JCC_CVaR.qY * system_info.xi;
        Joint_CVaR_Fy{j,i} = DRO_JCC_CVaR.fy;
        Joint_CVaR_FY{j,i} = DRO_JCC_CVaR.fY * system_info.xi;
            
        % Loop for each out-of-sample realization 
        for k = 1:OOS_sim
            system_info.Wreal = WPr(k,:)';
            system_info.DWreal = system_info.Wreal - system_info.mu;

            % Solve real-time optimal power flow for the solution of Bonferroni
            % approximation
            RT_solution_CVaR = RT_solve_R(system_info,DRO_ICC_CVaR.p,DRO_ICC_CVaR.ru,DRO_ICC_CVaR.rd);
            CVaR_Obj_IR(j,k,i) = RT_solution_CVaR.Obj_RT;  
            CVaR_lshed{j,k,i} = RT_solution_CVaR.lshed_RT;
            CVaR_flow{j,k,i} = DRO_ICC_CVaR.fy + DRO_ICC_CVaR.fY * system_info.DWreal;
            CVaR_p{j,k,i} = DRO_ICC_CVaR.Y * system_info.DWreal;
            CVaR_q{j,k,i} = DRO_ICC_CVaR.q + DRO_ICC_CVaR.qY * system_info.DWreal;
            CVaR_flag(j,k,i) = RT_solution_CVaR.Flag;

            
            % Solve real-time optimal power flow for the solution of Zymler
            % approximation
            RT_solution_Joint_CVaR = RT_solve_R(system_info,DRO_JCC_CVaR.p,DRO_JCC_CVaR.ru,DRO_JCC_CVaR.rd);
            Joint_CVaR_Obj_IR(j,k,i) = RT_solution_Joint_CVaR.Obj_RT;  
            Joint_CVaR_lshed{j,k,i} = RT_solution_Joint_CVaR.lshed_RT;
            Joint_CVaR_flow{j,k,i} = DRO_JCC_CVaR.fy + DRO_JCC_CVaR.fY * system_info.DWreal;
            Joint_CVaR_p{j,k,i} = DRO_JCC_CVaR.Y * system_info.DWreal;     
            Joint_CVaR_q{j,k,i} = DRO_JCC_CVaR.q + DRO_JCC_CVaR.qY * system_info.DWreal;
            Joint_CVaR_flag(j,k,i) = RT_solution_Joint_CVaR.Flag;

            
        end
        
        % Calculation of expected cost
        
        if DRO_ICC_CVaR.Flag == 0
            ICC_TC(j,i) = system_info.Cr1'*DRO_ICC_CVaR.ru + system_info.Cr2'*DRO_ICC_CVaR.rd + mean(CVaR_Obj_IR(j,:,i));
        else
            ICC_TC(j,i) = NaN;
        end
        
        if DRO_JCC_CVaR.Flag == 0
            JCC_TC(j,i) = system_info.Cr1'*DRO_JCC_CVaR.ru + system_info.Cr2'*DRO_JCC_CVaR.rd + mean(Joint_CVaR_Obj_IR(j,:,i));
        else
            JCC_TC(j,i) = NaN;
        end

    end
end

Time = toc
