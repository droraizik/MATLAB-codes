clc;
close all;
clear;

[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();


b1=backward_TM_abs_orig(:,10);
b2=backward_TM_abs_orig(:,20);


function score=PR_score1(b_col)
b_col=b_col(:)
score_num=mean(b_col.^2)-mean(b_col)^2
score_denom=(mean(b_col.^2)+1e2)
score=score_num/score_denom
fprintf('WWWWWWWWWWWWWWWWWWWWWWWWWWWWW\n');
end
