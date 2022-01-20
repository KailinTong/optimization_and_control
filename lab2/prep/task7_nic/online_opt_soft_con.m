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
global roh;
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
epsilon = sdpvar(1,1);
J = opt_var'*(Hy'*Qy*Hy + R)*opt_var + 2*opt_var'*Hy'*Qy*en + roh*epsilon;
global L;
w_bar = [L*un-u_min;u_max-L*un;gxn-x_min;x_max-gxn];
v_bar = ones(size(w_bar)); % all variables are soft constrained
constraints = [W*opt_var <= w_bar + v_bar*epsilon, 0 <= epsilon];

options = sdpsettings('verbose',0,'solver','qpoases','cachesolvers',1);
opt_diag = optimize(constraints,J,options);
if opt_diag.problem > 0
    error('Error during optimization')
else
    du = [value(opt_var(1));value(opt_var(2))];
end
end