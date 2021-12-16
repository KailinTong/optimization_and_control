clear all;
close all;
load('linearized_system.mat');

Q = [[10,0,0,0];
     [0,100,0,0];
     [0,0,0.01,0];
     [0,0,0,0.01];];
l = 0.3329;
mw = 2.9888;
V = 2.2696;
k1 = 11.1891;
ms = 0.5;
g = 9.81;

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    0, -ms*g/mw, -k1/mw, 0;
    0, (ms+mw)*g/(l*mw), k1/(l*mw), 0];
r = 0.1;
[kT, ~, E] = lqr(A, b, Q, r); 

