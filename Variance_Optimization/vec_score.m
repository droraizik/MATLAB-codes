function score=vec_score(b_col)
b_col=sqrt(b_col(:));
score_to_eps=mean(b_col.^2)-mean(b_col)^2;
score=score_to_eps/(mean(b_col.^2)+1e5);


% score=(mean(b_col.^2)-mean(b_col)^2)^2/mean(b_col);

