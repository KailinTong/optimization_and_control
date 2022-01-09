function [F,G,H,K] = calc_mpc_mat(A,B,C,Np,Nc,Q,R)
%% extract dimensions
n = size(A);
n = n(1);
m = size(B);
m = m(2);
p = size(C);
p = p(1);

%% calculate F
F = [];
for i=1:Np
    F = [F;C*A^i];
end

%% calculate G
G = [];
for i = 1:Np
    tempInnerMatrix = zeros(size(A));
    for j = 0:i-1
        tempInnerMatrix = tempInnerMatrix + A^j;
    end
    G = [G;C*tempInnerMatrix*B];
end

%% calculate H
H = [];
tempG = G;
zeroFiller = zeros(size(C*B));
for i=1:Nc
    H = [H, tempG];
    tempG = [zeroFiller; tempG];
    tempG = tempG(1:size(G,1),1:size(G,2));
end
    
%% calculate K
K = inv(H'*Q*H+R)*H'*Q;
selectMInputs = eye(m);
for i=1:(size(K,1)/m)-1
    selectMInputs = [selectMInputs, zeros(m)];
end
K = selectMInputs*K;

end

