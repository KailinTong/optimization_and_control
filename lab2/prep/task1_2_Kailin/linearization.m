syms x1 x2 u1 k1 k2 alpha1 beta1 gamma1
f = [-k1 * sqrt(x1) + alpha1 + sqrt(beta1 + gamma1 * u1), -k2 * sqrt(x2) + k1 * sqrt(x1)];

A = jacobian(f, [x1, x2])
b = jacobian(f, u1)

% A =
%  
% [ -k1/(2*x1^(1/2)),                0]
% [  k1/(2*x1^(1/2)), -k2/(2*x2^(1/2))]
%  
%  
% b =
%  
% gamma1/(2*(beta1 + gamma1*u1)^(1/2))
%                                     0

clear;

load('parameters.mat') 
load('equibrium.mat')
A = [-k1/(2*xe(1)^(1/2)), 0 ;
     k1/(2*xe(1)^(1/2)), -k2/(2*xe(2)^(1/2))];
disp(A)
b = [gamma1/(2*(beta1 + gamma1*ue(1))^(1/2)); 0];
disp(b)
save('linearized_matrix.mat', 'A', 'b')