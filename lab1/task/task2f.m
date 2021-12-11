clear all;
close all;

%% load measruements 4 and read data
load('measurements_4.mat');
t = measurements_4.time;
i_a = measurements_4.signals.values(:,1);
x_w = measurements_4.signals.values(:,2);
phi = measurements_4.signals.values(:,3); 
v_w = measurements_4.signals.values(:,4); 
omega = measurements_4.signals.values(:,5); 
omega_dot = measurements_4.signals.values(:,6); 

simin = [t, i_a];

load('L1_const.mat');
load('L2_const.mat');
load('Linf_const.mat');

%% L1: simulate
% set inital state and model constants
g = 9.81;
l = 30;
m_s = 0.5;
m_w = l1.m_w;
k_1 = l1.k_1;
V = l1.V;

sim('task2e_model_2015b');

figure(1);
subplot(2,2,1);
plot(t, x_w, 'DisplayName', 'data x_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,1), 'DisplayName', 'sim x_w');
xlabel('time');
ylabel('x_w');
title('L1 x_w');

subplot(2,2,2);
plot(t, phi, 'DisplayName', 'data \phi');
hold on;
plot(x_sim.time, x_sim.signals.values(:,2), 'DisplayName', 'sim \phi');
xlabel('time');
ylabel('\phi');
title('L1 \phi');

subplot(2,2,3);
plot(t, v_w, 'DisplayName', 'data v_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,3), 'DisplayName', 'sim v_w');
xlabel('time');
ylabel('v_w');
title('L1 v_w');

subplot(2,2,4);
plot(t, omega, 'DisplayName', 'data \omega');
hold on;
plot(x_sim.time, x_sim.signals.values(:,4), 'DisplayName', 'sim \omega');
xlabel('time');
ylabel('\omega');
title('L1 \omega');

%% L2: simulate
% set inital state and model constants
g = 9.81;
l = 30;
m_s = 0.5;
m_w = l2.m_w;
k_1 = l2.k_1;
V = l2.V;

sim('task2e_model_2015b');

figure(2);
subplot(2,2,1);
plot(t, x_w, 'DisplayName', 'data x_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,1), 'DisplayName', 'sim x_w');
xlabel('time');
ylabel('x_w');
title('L2 x_w');

subplot(2,2,2);
plot(t, phi, 'DisplayName', 'data \phi');
hold on;
plot(x_sim.time, x_sim.signals.values(:,2), 'DisplayName', 'sim \phi');
xlabel('time');
ylabel('\phi');
title('L2 \phi');

subplot(2,2,3);
plot(t, v_w, 'DisplayName', 'data v_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,3), 'DisplayName', 'sim v_w');
xlabel('time');
ylabel('v_w');
title('L2 v_w');

subplot(2,2,4);
plot(t, omega, 'DisplayName', 'data \omega');
hold on;
plot(x_sim.time, x_sim.signals.values(:,4), 'DisplayName', 'sim \omega');
xlabel('time');
ylabel('\omega');
title('L2 \omega');


%% Linf: simulate
% set inital state and model constants
g = 9.81;
l = 30;
m_s = 0.5;
m_w = linf.m_w;
k_1 = linf.k_1;
V = linf.V;

sim('task2e_model_2015b');

figure(3)
subplot(2,2,1);
plot(t, x_w, 'DisplayName', 'data x_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,1), 'DisplayName', 'sim x_w');
xlabel('time');
ylabel('x_w');
title('L infinity x_w');

subplot(2,2,2);
plot(t, phi, 'DisplayName', 'data \phi');
hold on;
plot(x_sim.time, x_sim.signals.values(:,2), 'DisplayName', 'sim \phi');
xlabel('time');
ylabel('\phi');
title('L infinity \phi');

subplot(2,2,3);
plot(t, v_w, 'DisplayName', 'data v_w');
hold on;
plot(x_sim.time, x_sim.signals.values(:,3), 'DisplayName', 'sim v_w');
xlabel('time');
ylabel('v_w');
title('L infinity v_w');

subplot(2,2,4);
plot(t, omega, 'DisplayName', 'data \omega');
hold on;
plot(x_sim.time, x_sim.signals.values(:,4), 'DisplayName', 'sim \omega');
xlabel('time');
ylabel('\omega');
title('L infinity \omega');
