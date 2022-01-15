load('linearized_matrix.mat')
sys = ss(A,b,[],[]);
Ts = 0.5;
sysd = c2d(sys,Ts);
disp(sysd.A);
disp(sysd.B)