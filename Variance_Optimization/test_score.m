clc;
close all;
clear;
[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();
num_elems=2;
element_location = create_element_matrix(N2,num_elems,[14,15]);
elems=find(element_location>0);
backward_TM=backward_TM_orig(:,elems);                                   %
backward_TM_abs=backward_TM_abs_orig(:,elems);                           %


at=[50e3;0];
vec_score(backward_TM_abs*at)
figure;imagesc(reshape(backward_TM_abs*at,N1,N1))


at=[40e3;0];
vec_score(backward_TM_abs*at)
figure;imagesc(reshape(backward_TM_abs*at,N1,N1))