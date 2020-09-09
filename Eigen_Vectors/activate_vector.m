function bt = activate_vector(bt,forward_TM,backward_TM_abs,forward_TM_show,elems_c,elems_r,N1,N2)
a_show=forward_TM_show*bt;                  %show a plane
figure;imagesc(abs(reshape(a_show,N2,N2)));
hold on;plot(elems_c,elems_r,'*r')
bt=bt/norm(bt);                             %normalization
at=abs(forward_TM*bt).^2;                   %intensity at beads
at=at/norm(at);                             %normalize bead intensity
bt=backward_TM_abs*at;                      %plane b in camera
