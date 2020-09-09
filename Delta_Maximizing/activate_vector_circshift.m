function camera_int = activate_vector_circshift(forward_TM,backward_TM,SLM_patt,N1,show_flag,forward_TM_show,N2,elems_c,elems_r)
num_elems=size(backward_TM,2);
flag = validate_SLM_patt(SLM_patt);
if(flag==1)
    warning('SLM abs is not zeros and ones');
end
%%
b_in=ifftshift(ifft2(ifftshift(SLM_patt)));             %vector to enter tissue
b_in_col=reshape(b_in,N1^2,1);                          %
b_in_col=b_in_col/norm(b_in_col);

if(show_flag==1)
    a_plane=forward_TM_show*b_in_col;
    a_plane=reshape(a_plane,N2,N2);
    figure;imagesc(abs(a_plane)); hold on
    plot(elems_c,elems_r,'*r');
end

a=abs(forward_TM*b_in_col).^2;                          %plane a
a=a/norm(a);
a2=repelem(a',N1^2,1);
b_out_cols=a2.*backward_TM;
b_out_images=reshape(b_out_cols,N1,N1,num_elems);
% figure;imagesc(abs(b_out_images))
b_out_images_F=fftshift(fft2(fftshift(b_out_images)));
final_field_F=conj(SLM_patt).*(b_out_images_F);      
final_field=ifftshift(ifft2(ifftshift(final_field_F))); %field at camera
camera_int=abs(final_field).^2; 
camera_int=sum(camera_int,3);
% figure;imagesc(sqrt(camera_int))




