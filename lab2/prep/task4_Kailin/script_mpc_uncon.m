clear all;
close all;

%% system constants
global k_1;
k_1 = 0.391;
global h_01;
h_01 = 10.17;
global alpha_1;
alpha_1 = 0.055;
global beta_1;
beta_1 = -3.077;
global gamma_1;
gamma_1 = 3.551;

global k_2;
k_2 = 0.386;
global h_02;
h_02 = 10.53;
global alpha_2;
alpha_2 = 0.059;
global beta_2; 
beta_2 = -3.246;
global gamma_2; 
gamma_2 = 3.610;

global u_max;
u_max = 5;
global u_min;
u_min = 0;

%% equilibrium
global ye;
ye = 15;
global xe2;
xe2 = 25.53;
global ue2;
ue2 = 0; 
global ue1;
ue1 = ((k_2*sqrt(xe2) - alpha_1)^2 - beta_1)/gamma_1; 
global xe1;
xe1 = (k_2*sqrt(xe2)/k_1)^2;
global xe;
xe = [xe1; xe2];
global ue;
ue = [ue1; ue2];

%% load models
load('../task2_nic/lin_ss_cont.mat');
load('../task3_nic/lin_ss_disc.mat');

%% model predictive control 
Np = 10;
Nc = 1;
Q = 1*eye(1*Np);
R = 1*eye(2*Nc);
[F, G, H, K] = cal_mpc_matrix(Ad, Bd, Cd, Np, Nc, Q, R);
[F_, G_, H_, K_] = calc_mpc_mat(Ad, Bd, Cd, Np, Nc, Q, R);

disp(sum(sum(F - F_)))
disp(sum(sum(G - G_)))
disp(sum(sum(H - H_)))
disp(sum(sum(K - K_)))
disp('calculate mpc matrix for both scripts is correct!')

%% simulation
time_end = 150;
r = [linspace(0,time_end*0.6,100*time_end*0.6)', repmat(repmat([ye*0.3],1,Np)',1,100*time_end*0.6)'];
r = [r;[linspace(time_end*0.6,time_end,100*time_end*0.4)', repmat(repmat([0.5*ye],1,Np)',1,100*time_end*0.4)']];

sim('mpc_uncon_Kailin');





