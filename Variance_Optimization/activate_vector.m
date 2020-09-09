function [b_out,at] = activate_vector(forward_TM,backward_TM_abs,b_in,show_flag,forward_TM_show,N2)
b_in=b_in(:);
b_in=b_in/norm(b_in);                             %normalization

if(show_flag==2)
    at=abs(forward_TM_show*b_in).^2;
    at=reshape(at,N2,N2);
    figure(123);imagesc(at);pause(0.001)

end

at=abs(forward_TM*b_in).^2;                   %intensity at beads
b_out=backward_TM_abs*at;                      %plane b in camera