clc;
clear;
close all;
%% Initialize
[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();
num_elems=10;

element_location_cell{1} = create_element_matrix(N2,num_elems,[10 12 14 16 28 30 32 34]);
element_location_cell{2} = create_element_matrix(N2,num_elems,[46 47 48 49 50 51 52 53]);
element_location_cell{3} = create_element_matrix(N2,num_elems,[10 19 28 29 30]);
element_location_cell{4} = create_element_matrix(N2,num_elems,[10 18 73 81]);
for elems_i=1:4
    for kk=1:3
        element_location = element_location_cell{elems_i};
        elems=find(element_location>0);
        num_elems=numel(elems);
        O_support=double(get_circular_mask(N1,N1,N1/8));                    %support for fourier field
        [elems_r,elems_c]=ind2sub([N2 N2],elems);                           %elemnt x,y
        forward_TM=forward_TM_show(elems,:);                                     %
        backward_TM=backward_TM_orig(:,elems);                                   %
        backward_TM_abs=backward_TM_abs_orig(:,elems);                           %
        best_H_SLM = optimize_SLM(forward_TM,N1,backward_TM_abs,O_support,forward_TM_show,N2);
        best_b_in=fftshift(fft2(fftshift(best_H_SLM)));
        best_b_in=best_b_in(:)/norm(best_b_in(:));
        a=forward_TM_show*(best_b_in);
        a=reshape(a,N2,N2);
        figure;imagesc(abs(a).^2);title('intensity at beads');hold on;colorbar
        plot(elems_c,elems_r,'*r');
    end
end