function du = online_opt(input)
du = [0;0];
%% get global vars
p = 1; % number of outputs
m = 2; % number of inputs
n = 2; % number of states
global Np;
global Nc;
global Hy;
global Qy;
global R;
global W;
global u_min;
global u_max;
global x_min;
global x_max;

%% read input
rn = input(1:p*Np);
gyn = input(p*Np+1:p*Np*2);
gxn = input(p*Np*2+1:p*Np*4);
un = input(p*Np*4+1:end);
en = gyn - rn;

%% optimization
yalmip('clear') 
opt_var = sdpvar(m*Nc,1); 
J = opt_var'*(Hy'*Qy*Hy + R)*opt_var + 2*opt_var'*Hy'*Qy*en;
global L;
global itcounter;
% to compensate for weird startup behaviour
% where first prediction can of the states can not be used
w_bar = [L*un-u_min;u_max-L*un;gxn-x_min;x_max-gxn];
constraints = [W*opt_var <= w_bar];
itcounter = itcounter + 1;

options = sdpsettings('verbose',0,'solver','qpoases','cachesolvers',1);
opt_diag = optimize(constraints,J,options);
if opt_diag.problem > 0
    error('Error during optimization')
else
    du = [value(opt_var(1));value(opt_var(2))];
end
end