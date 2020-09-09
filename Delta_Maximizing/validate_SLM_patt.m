function flag = validate_SLM_patt(SLM_patt)
ab=abs(SLM_patt);
check_zeros_ones=(ab>1e-5).*(abs(ab-1)>1e-5);
if(any(any(check_zeros_ones)))
    flag=1;
else
    flag=0;
end