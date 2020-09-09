function variab=normalize_var(variab)
var_col=reshape(variab,numel(variab),1);
var_col=var_col/norm(var_col);
variab=reshape(var_col,size(variab));