% Input data file

% Wind data

load('AV_AEMO');

% Scaling the system
Scale_Factor = 1;

% Electricity Network Data
% From       To     X       CAP 
ElNetwork=[
    101     102   0.0146    175; %l1
    101     103   0.2253    175; %l2
    101     105   0.0907    400; %l3
    102     104   0.1356    175; %l4
    102     106    0.205    175; %l5
    103     109   0.1271    400; %l6
    103     124    0.084    200; %l7
    104     109    0.111    175; %l8
    105     110    0.094    400; %l9
    106     110   0.0642    400; %l10
    107     108   0.0652    600; %l11
    108     109   0.1762    175; %l12
    108     110   0.1762    175; %l13
    109     111    0.084    200; %l14
    109     112    0.084    200; %l15
    110     111    0.084    200; %l16
    110     112    0.084    200; %l17
    111     113   0.0488    500; %l18
    111     114   0.0426    500; %l19
    112     113   0.0488    500; %l20
    112     123   0.0985    500; %l21
    113     123   0.0884    500; %l22
    114     116   0.0594    1000; %l23
    115     116   0.0172    500; %l24
    115     121   0.0249    1000; %l25
    115     124   0.0529    500; %l26
    116     117   0.0263    500; %l27
    116     119   0.0234    500; %l28
    117     118   0.0143    500; %l29
    117     122   0.1069    500; %l30
    118     121   0.0132    1000; %l31
    119     120   0.0203    1000; %l32
    120     123   0.0112    1000; %l33
    121     122   0.0692    500; %l34
    ];

% for l=1:34
%     if ElNetwork(l,4) <= 200
%         ElNetwork(l,4) = 250;
%     end
% end


% Calculation of PTDF matrix
PTDF_calc;

% Gen Data
% Pmin   Pmax    R+      R-     El_Bus    Cost      phi   pipe    
% GenDATA=[
%     0    152    40      40        1        12       12.65   1; %1
%     0    152    40      40        2        13       0       0; %2
%     0	 300    70      70        7        11       0       0; %3
%     0	 591    60      60        13       17       0       0; %4
%     0	  60    30      30        15       18       11.12   2; %5
%     0	 155    30      30        15       14       0       0; %6
%     0	 155    30      30        16       15       14.88   1; %7
%     0	 400    50      50        18       5        0       0; %8
%     0	 400    50      50        21       7        0       0; %9
%     0	 300    50      50        22       20       0       0; %10
%     0	 310    60      60        23       10.52    16.80   2; %11
%     0	 350    40      40        23       10.89    0       0; %12
%     ];

GenDATA=[
    0    152    40      40        1        12       12.65   1; %1
    0    152    40      40        2        13       13.45   3; %2
    0	 300    70      70        7        11       0       0; %3
    0	 591    60      60        13       17       0       0; %4
    0	  60    30      30        15       18       11.12   2; %5
    0	 155    30      30        15       14       0       0; %6
    0	 155    30      30        16       15       14.88   1; %7
    0	 400    50      50        18       5        0       0; %8
    0	 400    50      50        21       7        0       0; %9
    0	 300    50      50        22       20       0       0; %10
    0	 310    60      60        23       10.52    16.80   2; %11
    0	 350    40      40        23       10.89    15.60   3; %12
    ];

% Pipeline capacity 
PipeCap = [10000; 5500; 7000];

% Demand       El_Bus  Cshed           
DemandDATA=[
    0.038      1     2000; %d1
    0.034      2     2000; %d2
    0.063      3     2000; %d3
    0.026      4     2000; %d4
    0.025      5     2000; %d5
    0.048      6     2000; %d6
    0.044      7     2000; %d7
     0.06      8     2000; %d8
    0.061      9     2000; %d9
    0.068      10    2000; %d10
    0.093      13    2000; %d11
    0.068      14    2000; %d12
    0.111      15    2000; %d13
    0.035      16    2000; %d14
    0.117      18    2000; %d15
    0.064      19    2000; %d16
    0.045      20    2000; %d17
    ];

DemandDATA(:,3) = 1000;

Total_Demand = 2650;

% Wmin   Wmax   El_Bus
WindDATA=[
    0     200     1; %1 
    0     200     2; %2 
    0     200     3; %3 
    0     200     4; %4 
    0     200     5; %5 
    0     200     6; %6 
%     0     200     7; %7 
%     0     200     8; %8
%     0     200     9; %1 
%     0     200     10; %2 
%     0     200     11; %3 
%     0     200     12; %4 
%     0     200     13; %5 
%     0     200     14; %6 
%     0     200     15; %7 
%     0     200     16; %8 
%     0     200     17; %3 
%     0     200     18; %4 
%     0     200     19; %5 
%     0     200     20; %6 
%     0     200     21; %7 
%     0     200     22; %8
    ];

WindDATA(:,2) = 250;
WindDATA(:,3) = [1,2,11,12,12,16];
% WindDATA(:,3) = [1,2,11,12,11,12,16,20];

% Fill system info structure
system_info = [];

system_info.PTDF = round([PTDF_nrf(:,1:ref_node-1), zeros(size(PTDF_nrf,1),1), PTDF_nrf(:,ref_node:size(PTDF_nrf,2))],2); % The final PTDF matrix
system_info.F = ElNetwork(:,4)/Scale_Factor;
system_info.D = Total_Demand * DemandDATA(:,1)/Scale_Factor;
system_info.Pmax = GenDATA(:,2)/Scale_Factor;
system_info.Pmin = GenDATA(:,1)/Scale_Factor;
system_info.R = (system_info.Pmax + system_info.Pmin)/2;
system_info.ResCap = GenDATA(:,2)*0.40;

system_info.FP = PipeCap;

system_info.Wmax = WindDATA(:,2)/Scale_Factor;
system_info.DiagWmax = diag(system_info.Wmax);
system_info.Wmin = WindDATA(:,1)/Scale_Factor;

%system_info.C = 0.5*round(GenDATA(:,6)*Scale_Factor,3);
system_info.Clsh = DemandDATA(:,3)*Scale_Factor;
system_info.C = [35;40;30;55;60;45;50;10;15;65;20;25]/2;
% system_info.Cr = [60;70;50;100;110;80;90;10;20;120;30;40];
system_info.Cr1 = 0.2*system_info.C;
system_info.Cr2 = 0.2*system_info.C;

% Mapping on the network
system_info.AG = zeros(N_El_nodes, size(GenDATA,1));
system_info.AW = zeros(N_El_nodes, size(WindDATA,1));
system_info.AD = zeros(N_El_nodes, size(DemandDATA,1));
system_info.PG = zeros(size(PipeCap,2), size(GenDATA,1));

for n=1:N_El_nodes
    for gg=1:size(GenDATA,1)
        if GenDATA(gg,5) == n
        system_info.AG(n,gg) = 1;
        end
    end
    for ww=1:size(WindDATA,1)
        if WindDATA(ww,3) == n
        system_info.AW(n,ww) = 1;
        end
    end
    for dd=1:size(DemandDATA,1)
        if DemandDATA(dd,2) == n
        system_info.AD(n,dd) = 1;
        end
    end
end

for pp = 1:size(PipeCap,1)
    for ggg = 1:size(GenDATA,1)
        if GenDATA(ggg,8) == pp
            system_info.PG(pp,ggg) = 1;
        end
    end
end

system_info.PG = system_info.PG .* GenDATA(:,7)';

% PTDF and mapping matrices combination

system_info.Qg = system_info.PTDF*system_info.AG;
system_info.Qw = system_info.PTDF*system_info.AW;
system_info.Qd = system_info.PTDF*system_info.AD;

