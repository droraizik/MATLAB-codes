function show_a_plane(forward_TM_show,SLM_pattern,N1,N2,elems_c,elems_r,plot_title)
bt=ifftshift(ifft2(ifftshift(SLM_pattern)));
bt=reshape(bt,N1^2,1);
bt=bt/norm(bt);
at=forward_TM_show*bt;
at=reshape(at,N2,N2);
figure;imagesc(abs(at));title(plot_title);
hold on
plot(elems_c,elems_r,'*r');