function dx = nonlinear_2_tank_model_function(input)
dx = [0;0;];

u=input(1:2);
x=input(3:4);

global k_1;
global alpha_1;
global beta_1;
global gamma_1;

global k_2;
global alpha_2;
global beta_2; 
global gamma_2; 

dx(1) = -k_1*sqrt(x(1)) + z_pi(alpha_1, beta_1, gamma_1, u(1));
dx(2) = -k_2*sqrt(x(2)) + k_1*sqrt(x(1)) + z_pi(alpha_2, beta_2, gamma_2, u(2));

end
