function score=PR_score(b_col)
b_col=b_col(:);
score=sqrt(mean(b_col.^2)-mean(b_col)^2);
score=score/(sqrt(mean(b_col.^2)+1e5));

