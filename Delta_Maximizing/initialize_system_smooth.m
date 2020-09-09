function [orig_f_TM,abs_b_TM,orig_b_TM,N1,N2]=initialize_system_smooth()
load('C:\Users\droraizik\Dropbox\Matrices\TM_test.mat','T');
u=T;
N1=size(u,1);
N2=size(u,3);


abs_b_TM=abs(reshape(u,N1*N1,N2*N2)).^2;                    % backward TM - abs value
element_strength_factor=50000/(mean(sum(abs_b_TM)));  
abs_b_TM=element_strength_factor*abs_b_TM;

orig_b_TM=(reshape(u,N1*N1,N2*N2));                         %backward TM - original value
element_strength_factor=50000/(mean(sum(orig_b_TM)));  
orig_b_TM=element_strength_factor*orig_b_TM;

orig_f_TM=orig_b_TM.';
% orig_f_TM=orig_b_TM';