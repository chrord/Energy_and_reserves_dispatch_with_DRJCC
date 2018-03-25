function [ CC_m ] = CC_matrices(si, DRO_param)

% Matrices for generation 

A_g = [zeros(size(si.Pmax,1),size(si.Pmax,1)),-eye(size(si.Pmax,1),size(si.Pmax,1)),zeros(size(si.Pmax,1),size(si.Pmax,1));
       zeros(size(si.Pmax,1),size(si.Pmax,1)),zeros(size(si.Pmax,1),size(si.Pmax,1)),-eye(size(si.Pmax,1),size(si.Pmax,1))];
   
B_g = [eye(size(si.Pmax,1),size(si.Pmax,1));
       -eye(size(si.Pmax,1),size(si.Pmax,1))];  
   
C_g = [zeros(size(si.Pmax,1),size(si.Wmax,1));
       zeros(size(si.Pmax,1),size(si.Wmax,1))]; 
   
d_g = [zeros(size(si.Pmax,1),1);
       zeros(size(si.Pmax,1),1)];   

% Matrices for transmission lines 

A_l = [si.Qg,zeros(size(si.F,1),size(si.Pmax,1)),zeros(size(si.F,1),size(si.Pmax,1));
       -si.Qg,zeros(size(si.F,1),size(si.Pmax,1)),zeros(size(si.F,1),size(si.Pmax,1))];
   
B_l = [si.Qg;
       -si.Qg];  
   
C_l = [si.Qw*si.DiagWmax;
       -si.Qw*si.DiagWmax]; 
   
d_l = [si.F - si.Qw*si.DiagWmax*si.mu + si.Qd*si.D;
       si.F + si.Qw*si.DiagWmax*si.mu - si.Qd*si.D];  
   
% Matrices for pipelines 

A_m = [si.PG,zeros(size(si.FP,1),size(si.Pmax,1)),zeros(size(si.FP,1),size(si.Pmax,1));
       -si.PG,zeros(size(si.FP,1),size(si.Pmax,1)),zeros(size(si.FP,1),size(si.Pmax,1))];
   
B_m = [si.PG;
       -si.PG];  
   
C_m = [zeros(size(si.FP,1),size(si.Wmax,1));
       zeros(size(si.FP,1),size(si.Wmax,1))]; 
   
d_m = [si.FP;
       zeros(size(si.FP,1),1)];  
   
   
jcc{1,1} = A_g;
jcc{2,1} = A_l;
jcc{3,1} = A_m;

jcc{1,2} = B_g;
jcc{2,2} = B_l;
jcc{3,2} = B_m;

jcc{1,3} = C_g;
jcc{2,3} = C_l;
jcc{3,3} = C_m;

jcc{1,4} = d_g;
jcc{2,4} = d_l;
jcc{3,4} = d_m;

jcc{1,5} = DRO_param.eps_joint_cvar;
jcc{2,5} = DRO_param.eps_joint_cvar;
jcc{3,5} = DRO_param.eps_joint_cvar;

CC_m.jcc = jcc;

end

