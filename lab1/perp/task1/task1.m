clear all;
close all;
clc;

%% Load and plot recorded data
load('measurements_1.mat');
% x1 = phi
% x2 = w
% x3 = dw/dt


% dx1/dt = x2, x3 = g/l * sin(x1)
t = measurements_1.time;
x1 = measurements_1.signals.values(:,1);
x2 = measurements_1.signals.values(:,2);
x3 = measurements_1.signals.values(:,3); 
figure() 
subplot(2,2,1)
plot(t, x1, '-b', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
subplot(2,2,2)
plot(t, x2, '-b', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
subplot(2,2,3)
plot(t, x3 , '-b', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_3'), grid on;

%% identification (L1 - optimization)
p = sdpvar(1,1); % p = g/l
e = x3 - p(1)*sin(x1);

bounds = sdpvar(length(e),1);
constraints = [-bounds <= e <= bounds, p(1) >= 0];
options = sdpsettings('verbose',0);
diagnostics_L1 = optimize(constraints,sum(bounds),options);
if diagnostics_L1.problem > 0
error('Error during optimization')
end
p_hat = value(p);
g = 9.8;
l_hat_l1 = g / p_hat(1);


%% identification (L2 - optimization)
p = sdpvar(1,1); % p = g/l
constraints = [p(1) >= 0];

e = x3 - p(1)*sin(x1);

options = sdpsettings('verbose', 0);
diagnostics_L2 = optimize(constraints, e'*e, options);
if diagnostics_L2.problem > 0
    error('Error during L2-Optimization')
end
p_hat = value(p);
g = 9.8;
l_hat_l2 = g / p_hat(1);


%% identification (L infinite - optimization)
p = sdpvar(1,1); % p = g/l
e = x3 - p(1)*sin(x1);

bounds = sdpvar(1,1);
constraints = [-bounds <= e <= bounds, p(1) >= 0];
options = sdpsettings('verbose',0);
diagnostics_Linf = optimize(constraints,bounds,options);
if diagnostics_Linf.problem > 0
error('Error during optimization')
end
p_hat = value(p);
g = 9.8;
l_hat_linf = g / p_hat(1);



Table_parameters = table([1;l_hat_l1], [2;l_hat_l2], [3;l_hat_linf], 'VariableNames', {'l_hat_by_L1', 'l_hat_by_L2', 'l_hat_by_Linf'}, 'RowNames', {'true', 'identified'});
disp(Table_parameters);



%% Run simulation and compare results
Tsim = t(end); Ts = 0.01;
x0 = [x1(1); x2(1)];
l = l_hat_l1;
sim('task1c');
% 
figure()
subplot(3,2,1), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 1', 'simulation L1');
subplot(3,2,2), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 1', 'simulation L1');

l = l_hat_l2;
sim('task1c');
subplot(3,2,3), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 1', 'simulation L2');
subplot(3,2,4), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 1', 'simulation L2');

l = l_hat_linf;
sim('task1c');
subplot(3,2,5), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 1', 'simulation Linf');
subplot(3,2,6), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 1', 'simulation Linf');

 
%% Load and  data
load('measurements_2.mat');
% x1 = phi
% x2 = w
% x3 = dw/dt


% dx1/dt = x2, x3 = g/l * sin(x1)
t = measurements_2.time;
x1 = measurements_2.signals.values(:,1);
x2 = measurements_2.signals.values(:,2);
x3 = measurements_2.signals.values(:,3); 
x0 = [x1(1); x2(1)];
l = l_hat_l1;
sim('task1c');

% 
figure()
subplot(3,2,1), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 2', 'simulation L1');
subplot(3,2,2), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 2', 'simulation L1');

l = l_hat_l2;
sim('task1c');
subplot(3,2,3), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 2', 'simulation L2');
subplot(3,2,4), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 2', 'simulation L2');

l = l_hat_linf;
sim('task1c');
subplot(3,2,5), hold on;
plot(t, x1, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,1), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_1'), grid on;
legend('recorded 2', 'simulation Linf');
subplot(3,2,6), hold on;
plot(t, x2, '-b', 'LineWidth', 1.5)
plot(x_sim.time, x_sim.signals.values(:,2), ':r', 'LineWidth', 1.5)
xlabel('t in s'), ylabel('x_2'), grid on;
legend('recorded 2', 'simulation Linf');





