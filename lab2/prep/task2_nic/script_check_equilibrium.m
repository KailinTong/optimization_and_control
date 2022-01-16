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

%% linearization

A = [[(-k_1/(2*sqrt(xe1))), 0]; [(k_1/(2*sqrt(xe1))), (-k_2/(2*sqrt(xe2)))]];
B = [[(gamma_1/(2*sqrt(beta_1+gamma_1*ue1))), 0];[0, 0]];
C = [0,1];
D = [0,0];

%% simulation
time_end = 10;
simin = [linspace(0,time_end,100*time_end)', repmat([ue1],1,100*time_end)', repmat([ue2],1,100*time_end)'];

sim('check_equilibrium');





