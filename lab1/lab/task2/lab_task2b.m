clear all;
close all;

%% Load and plot recorded data
load('lab_measurement_3.mat');

start_index = 150;
t = x2_meas.time(start_index:end);
i_a = u2_meas.signals.values(start_index:end,1);
x_w = x2_meas.signals.values(start_index:end,1);
phi = x2_meas.signals.values(start_index:end,2); 
v_w = x2_meas.signals.values(start_index:end,3); 
omega = x2_meas.signals.values(start_index:end,4); 
omega_dot = dw2_meas.signals.values(start_index:end); 

% known parameters
l = 0.3329;
g = 9.8;


y = omega_dot*l - g*sin(phi);
f1 = omega_dot.*(cos(phi)).^2 - omega.^2.*sin(phi).*cos(phi);
f2 = -cos(phi).*i_a;
f3 = cos(phi).*v_w;

%% identification (L1 - optimization)
yalmip('clear') % delete old yalmip variables
p = sdpvar(1,3); 
e = y - p(1)*f1 - p(2)*f2 - p(3)*f3;
bounds = sdpvar(length(e),1);
constraints = [-bounds <= e <= bounds, p(1)>=0, p(2)>=0, p(3)>=0];

options = sdpsettings('verbose',0);
diagnostics_L1 = optimize(constraints,sum(bounds),options);
if diagnostics_L1.problem > 0
error('Error during optimization')
end
p_l1 = value(p);
disp('L1 estimated parameters:')
disp(p_l1)


%% identification (L2 - optimization)
yalmip('clear') % delete old yalmip variables
p = sdpvar(1,3); 
e = y - p(1)*f1 - p(2)*f2 - p(3)*f3;
constraints = [p(1)>=0, p(2)>=0, p(3)>=0];

options = sdpsettings('verbose', 0);
diagnostics_L2 = optimize(constraints, e'*e, options);
if diagnostics_L2.problem > 0
    error('Error during L2-Optimization')
end
p_l2 = value(p);
disp('L2 estimated parameters:')
disp(p_l2)



%% identification (L infinite - optimization)
yalmip('clear') % delete old yalmip variables
p = sdpvar(1,3); 
e = y - p(1)*f1 - p(2)*f2 - p(3)*f3;
bounds = sdpvar(1,1);
constraints = [-bounds <= e <= bounds, p(1)>=0, p(2)>=0, p(3)>=0];

options = sdpsettings('verbose',0);
diagnostics_Linf = optimize(constraints,bounds,options);
if diagnostics_Linf.problem > 0
error('Error during optimization')
end
p_linf = value(p);
disp('L infinite  estimated parameters:')
disp(p_linf)


save('task2b_p.mat', 'p_l1', 'p_l2', 'p_linf');








