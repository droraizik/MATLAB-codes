clc;
close all;
clear;

[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();


b1=backward_TM_abs_orig(:,10);
b2=backward_TM_abs_orig(:,20);

b_a=b1+b2;
b_a_1s=b_a/norm(b_a);
b_b=b1;
b_b_1s=b_b/norm(b_b);

var(b_a)
var(b_a_1s)
var(b_b)
var(b_b_1s)




function score=PR_score1(b_col)
b_col=b_col(:)
score_num=mean(b_col.^2)-mean(b_col)^2
score_denom=(mean(b_col.^2)+1e2)
score=score_num/score_denom
fprintf('WWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n');
end
