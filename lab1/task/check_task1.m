clear all;
close all;

%% load data
load('measurements_1.mat');
t = measurements_1.time;
phi = measurements_1.signals.values(:,1);
omega = measurements_1.signals.values(:,2);
omega_dot = measurements_1.signals.values(:,3); 


%% definitions
g = 9.81;

%% check
% l = (g/omega_dot)*sin(phi)
l = (g./omega_dot).*sin(phi);

% plot result
figure(1);
plot(t,l);
title("task 1: length verification");
xlabel("time t");
ylabel("length l");

