clear; clc;

k1 = 0.391;
k2 = 0.386;
h01 = 10.17;
h02 = 10.53;
alpha1 = 0.055;
alpha2 = 0.059;
beta1 = -3.077;
beta2 = -3.246;
gamma1 = 3.551;
gamma2 = 3.610;
save('parameters.mat')

ye = 15;
he2 = ye;
xe2 = he2 + h02;
xe1 = (k2 / k1)^2 * xe2;
he1 = xe1 - h01;
ue1 = ((k1 * sqrt(xe1) - alpha1)^2 - beta1) / gamma1;
ue2 = 0;

xe = [xe1; xe2];
ue = [ue1; ue2];
save('equibrium.mat', 'xe', 'ue');


