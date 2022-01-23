function [F, G, H, K] = cal_mpc_matrix(A, B, C, Np, Nc, Q, R)
%CAL_MPC_MATRIX calcuate the matrix for the MPC unconstrained case
%   Detailed explanation goes here

n = size(A);
n = n(1);
m = size(B);
m = m(2);
p = size(C);
p = p(1);

F = [];
for i = 1:Np
    F = [F; C*A^i];
end

G = [];
for i = 1:Np
    innter_mat = eye(n);
    for j = 1:i-1
        innter_mat = innter_mat + A^j;
    end
    G = [G; C * innter_mat * B];
end

H = []
for i = 1:Nc
    column = zeros(size(G));
    column((i-1)*p+1:end,1:end) = G(1:p*(Np-i+1),1:end);
    H = [H, column];
end


K_full = inv((H'*Q*H + R))*H'*Q;
size_K_full = size(K_full);
zero_filter = zeros(m, size_K_full(1));
zero_filter(1:m, 1:m) = eye(m);
K = zero_filter * K_full;
end

