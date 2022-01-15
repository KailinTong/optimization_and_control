function dx = nonlinear_2_tank_model_function(input)
dx = [0;0;];

u=input(1:2);
x=input(3:4);

global k1;
global alpha1;
global beta1;
global gamma1;

global k2;
global alpha2;
global beta2; 
global gamma2; 

dx(1) = -k1*sqrt(x(1)) + z_pi(alpha1, beta1, gamma1, u(1));
dx(2) = -k2*sqrt(x(2)) + k1*sqrt(x(1)) + z_pi(alpha2, beta2, gamma2, u(2));

end
