%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
% Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Node data
N_El_nodes = 24;
ref_node = 13;

% ADMITANCE MATRIX
B_N = zeros(N_El_nodes, N_El_nodes);

% Off-diagonal elements B-matrix
for l=1:size(ElNetwork,1) % Number of Lines
    B_N(ElNetwork(l,1)-100, ElNetwork(l,2)-100) = -1/ElNetwork(l,3);
    B_N(ElNetwork(l,2)-100, ElNetwork(l,1)-100) = -1/ElNetwork(l,3);
end

% Diagonal elements B-matrix
for k=1:N_El_nodes
    B_N(k,k) = -sum(B_N(k,:));
end

B_L = zeros(size(ElNetwork,1), N_El_nodes);

% Off-diagonal elements B-matrix
for l=1:size(ElNetwork,1) % Number of Lines
    B_L(l, ElNetwork(l,1)-100) = 1/ElNetwork(l,3);
    B_L(l, ElNetwork(l,2)-100) = -1/ElNetwork(l,3);
end

% Remove ref node
B_NN = B_N;
B_NN(:,ref_node)=[];
B_NN(ref_node,:)=[];

B_LL = B_L;
B_LL(:,ref_node) = [];

PTDF_nrf = B_LL * (B_NN ^ (-1));