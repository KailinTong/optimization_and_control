function dx = linear_2_tank_model_function(input)
dx = [0;0;];

u=input(1:2);
x=input(3:4);

global Jx;
global Ju;

dx = Jx*x + Ju*u;

end
