input = [15;15;15;15;15;15;15;15;15;15;0.210346529891715;0.479304083419933;0.853104446335800;1.32763096691154;1.89888626819646;2.56298918827729;3.31617179401091;4.15477646653714;5.07525305691677;6.07415611027762;11.1831018684675;0.210346529891715;16.9023781233441;0.479304083419933;22.5106670270429;0.853104446335800;28.0101223821733;1.32763096691154;33.4028561949991;1.89888626819646;38.6909394865319;2.56298918827729;43.8764030878836;3.31617179401091;48.9612384201861;4.15477646653714;53.9473982593742;5.07525305691677;58.8367974861285;6.07415611027762;12.7986063403360;-5.55111512312578e-17];
du = [0;0];
disp('break');
%% get global vars
p = 1; % number of outputs
m = 2; % number of inputs
n = 2; % number of states
global Np;
global Nc;
global Hy;
global Qy;
global R;
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
w_bar = [L*un-u_min;u_max-L*un;gxn-x_min;x_max-gxn];
constraints = [W*opt_var <= w_bar];

options = sdpsettings('verbose',0,'solver','qpoases','cachesolvers',1);
opt_diag = optimize(constraints,J,options);
if opt_diag.problem > 0
    error('Error during optimization')
end
du = [value(opt_var(1));value(opt_var(2))];