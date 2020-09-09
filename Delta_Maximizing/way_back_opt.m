clc;
clear;
close all;
[forward_TM,backward_TM_abs,backward_TM,N1,N2]=initialize_system_smooth();
num_elems=2;
element_location = create_element_matrix(N2,num_elems,[37, 46,55,56,57]);
elems=find(element_location>0);
num_elems=numel(elems);
O_support=double(get_circular_mask(N1,N1,N1/8));                    %support for fourier field
[elems_r,elems_c]=ind2sub([N2 N2],elems);                           %elemnt x,y
forward_TM_show=forward_TM;                                         %
forward_TM=forward_TM(elems,:);                                     %
backward_TM=backward_TM(:,elems);                                   %
backward_TM_abs=backward_TM_abs(:,elems);                           %

a_plane=zeros(numel(elems),1);
a_plane(1)=1;
a_plane(2)=1;
a_plane(3)=1;
a_plane(4)=1;
a_plane(5)=1;

%% show a-plane
a_plane_show=zeros(N2^2,1);
a_plane_show(elems(:))=1;
a_plane_show=reshape(a_plane_show,N2,N2);
figure;imagesc(abs(a_plane_show)); hold on
plot(elems_c,elems_r,'*r');
%%

a_plane=a_plane/norm(a_plane);
a2=repelem(a_plane',N1^2,1);
b_out_cols=a2.*backward_TM;
b_out_images=reshape(b_out_cols,N1,N1,num_elems);
b_out_images_F=fftshift(fft2(fftshift(b_out_images)));
SLM_patt=conj(sign(b_out_images_F(:,:,1)));

final_field_F=(SLM_patt).*(b_out_images_F);      
final_field=ifftshift(ifft2(ifftshift(final_field_F))); %field at camera
camera_int=abs(final_field).^2; 
camera_int=sum(camera_int,3);
figure;imagesc(rot90(sqrt(camera_int),3))

