clear all;
close all;
load('linearized_system.mat');

Q = [[1,0,0,0];
     [0,100,0,0];
     [0,0,1,0];
     [0,0,0,100];];
r = 0.1;
[K, ~, E] = lqr(A, b, Q, r); 

