clear all;
close all;

load quadcopter_data_unc.mat
epsilon = 0.00001;

%% call calc function
[myF, myG, myH, myK] = cal_mpc_matrix(A,B,C,nP,nC,Q,R);

%% test F
if all(all(myF >= F-epsilon)) && all(all(myF <= F+epsilon))
    disp('F is correct');
else
    disp('F is false');
end

%% test G
if all(all(myG >= G-epsilon)) && all(all(myG <= G+epsilon))
    disp('G is correct');
else
    disp('G is false');
end

%% test H
if all(all(myH >= H-epsilon)) && all(all(myH <= H+epsilon))
    disp('H is correct');
else
    disp('H is false');
end

%% test K
if all(all(myK >= K-epsilon)) && all(all(myK <= K+epsilon))
    disp('K is correct');
else
    disp('K is false');
end