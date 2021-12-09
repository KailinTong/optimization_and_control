load('task2b_p.mat');

% known parameters
l = 0.3;
g = 9.8;
m_s = 0.5;


l1 = struct('m_w',0,'V',0,'k_1',0); 
l2 = struct('m_w',0,'V',0,'k_1',0); 
linf = struct('m_w',0,'V',0,'k_1',0);

l1.m_w = (l/p_l1(1) - 1)*m_s;
l1.V =  (l1.m_w + m_s)*p_l1(2); 
l1.k_1 = (l1.m_w + m_s)*p_l1(3); 
disp('L1  estimated parameters:')
disp(l1)
save('L1_const.mat', 'l1');

l2.m_w = (l/p_l2(1) - 1)*m_s;
l2.V =  (l2.m_w + m_s)*p_l2(2); 
l2.k_1 = (l2.m_w + m_s)*p_l2(3); 
disp('L2  estimated parameters:')
disp(l2)
save('L2_const.mat', 'l2');

linf.m_w = (l/p_linf(1) - 1)*m_s; 
linf.V =  (linf.m_w + m_s)*p_linf(2); 
linf.k_1 = (linf.m_w + m_s)*p_linf(3); 
disp('L infinite  estimated parameters:')
disp(linf)
save('Linf_const.mat', 'linf');
