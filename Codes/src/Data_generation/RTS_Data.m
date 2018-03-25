%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data file

% Wind data from ref. [31]
load('AV_AEMO');

% Electricity network data
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


% Calculation of PTDF matrix
PTDF_calc;

% Power production units data
GenDATA=[
    0    152        1       12.65   1; %1
    0    152        2       13.45   3; %2
    0	 300        7       0       0; %3
    0	 591        13      0       0; %4
    0	  60        15      11.12   2; %5
    0	 155        15      0       0; %6
    0	 155        16      14.88   1; %7
    0	 400        18      0       0; %8
    0	 400        21      0       0; %9
    0	 300        22      0       0; %10
    0	 310        23      16.80   2; %11
    0	 350        23      15.60   3; %12
    ];

% Pipeline capacity 
PipeCap = [10000; 5500; 7000];

% Electricity demand data
% Demand       El_Bus  Cshed           
DemandDATA=[
    0.038      1     1000; %d1
    0.034      2     1000; %d2
    0.063      3     1000; %d3
    0.026      4     1000; %d4
    0.025      5     1000; %d5
    0.048      6     1000; %d6
    0.044      7     1000; %d7
     0.06      8     1000; %d8
    0.061      9     1000; %d9
    0.068      10    1000; %d10
    0.093      13    1000; %d11
    0.068      14    1000; %d12
    0.111      15    1000; %d13
    0.035      16    1000; %d14
    0.117      18    1000; %d15
    0.064      19    1000; %d16
    0.045      20    1000; %d17
    ];

Total_Demand = 2650;

% Wind farm data
% Wmin   Wmax   El_Bus
WindDATA=[
    0     250     1; %1 
    0     250     2; %2 
    0     250     11; %3 
    0     250     12; %4 
    0     250     12; %5 
    0     250     16; %6 
    ];

% Fill system info structure
system_info = [];

system_info.PTDF = round([PTDF_nrf(:,1:ref_node-1), zeros(size(PTDF_nrf,1),1), PTDF_nrf(:,ref_node:size(PTDF_nrf,2))],2); % The final PTDF matrix
system_info.F = ElNetwork(:,4);
system_info.D = Total_Demand * DemandDATA(:,1);
system_info.Pmax = GenDATA(:,2);
system_info.Pmin = GenDATA(:,1);
system_info.R = (system_info.Pmax + system_info.Pmin)/2;
system_info.ResCap = GenDATA(:,2)*0.40;
system_info.FP = PipeCap;
system_info.Wmax = WindDATA(:,2);
system_info.DiagWmax = diag(system_info.Wmax);
system_info.Wmin = WindDATA(:,1);
system_info.Clsh = DemandDATA(:,3);
system_info.C = [35;40;30;55;60;45;50;10;15;65;20;25]/2;
system_info.Cr1 = 0.2*system_info.C;
system_info.Cr2 = 0.2*system_info.C;

% Mapping on the network
system_info.AG = zeros(N_El_nodes, size(GenDATA,1));
system_info.AW = zeros(N_El_nodes, size(WindDATA,1));
system_info.AD = zeros(N_El_nodes, size(DemandDATA,1));
system_info.PG = zeros(size(PipeCap,2), size(GenDATA,1));

for n=1:N_El_nodes
    for gg=1:size(GenDATA,1)
        if GenDATA(gg,3) == n
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
        if GenDATA(ggg,5) == pp
            system_info.PG(pp,ggg) = 1;
        end
    end
end

system_info.PG = system_info.PG .* GenDATA(:,4)';

% PTDF and mapping matrices combination
system_info.Qg = system_info.PTDF*system_info.AG;
system_info.Qw = system_info.PTDF*system_info.AW;
system_info.Qd = system_info.PTDF*system_info.AD;

