clc;
clear;
close all;
%% Initialize
num_elems=1;
[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();
element_location_cell{1} = create_element_matrix(N2,num_elems,[10 12 14 16 28 30 32 34]);
element_location_cell{2} = create_element_matrix(N2,num_elems,[46 47 48 49 50 51 52 53]);
element_location_cell{3} = create_element_matrix(N2,num_elems,[10 19 28 29 30]);
element_location_cell{4} = create_element_matrix(N2,num_elems,[10 18 73 81]);
for elems_i=1:1
    for reps_i=1:3
        element_location=element_location_cell{elems_i};
        elems=find(element_location>0);
        num_elems=numel(elems);O_support=double(get_circular_mask(N1,N1,N1/8));                    %support for fourier field
        [elems_r,elems_c]=ind2sub([N2 N2],elems);                           %elemnt x,y
        forward_TM=forward_TM_show(elems,:);                                     %
        backward_TM=backward_TM_orig(:,elems);                                   %
        backward_TM_abs=backward_TM_abs_orig(:,elems);                           %
        PR_tries=2;
        big_it=4;
        init_tries=20;
        delta_ints=zeros(1,PR_tries*2+1);
        SLM_patterns=zeros(N1,N1,PR_tries*2+1);
        init_patterns=zeros(N1,N1,init_tries);
        init_tries_deltas=zeros(1,init_tries);
        SLM_patts_init=zeros(N1,N1,PR_tries*2+1);
        [X,Y]=meshgrid(1:N1);
        sigma=4;
        z=exp(-0.5/sigma.*((X-N1/2).^2+(Y-N1/2).^2));
        fprintf("Number of Camera Shots: %d\n",init_tries*PR_tries*2*2+big_it*PR_tries*2);
        %% First try
        for init_i=1:init_tries
            fprintf("Initialization try: PR try %d\n",init_i);
            SLM_patt=conv2(exp(1i*2*pi*rand(N1)),z,'same');  %first try - smooth rand phase
            SLM_patt=sign(SLM_patt).*O_support;
            camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
            for i=1:PR_tries*2
%                 fprintf("Initialization SLM: PR try %d\n",i);
                %% Original solution
                camera_field_PR=GS_PR_max_mid(camera_int,N1);
                SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR))));
                SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                camera_int=normalize_var(camera_int);
                delta_ints(2*i-1)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                SLM_patts_init(:,:,i*2-1)=SLM_patt_new;
                %% Conjugate solution
                camera_field_PR_conj=conj(camera_field_PR);
                SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR_conj))));
                SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                camera_int=normalize_var(camera_int);
                delta_ints(2*i)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                SLM_patts_init(:,:,i*2)=SLM_patt_new;
            end
            [init_tries_deltas(init_i),max_ind]=max(delta_ints);
            init_patterns(:,:,init_i)=SLM_patts_init(:,:,max_ind);

        end
        [~,max_init_patt]=max(init_tries_deltas);
        SLM_patt=init_patterns(:,:,max_init_patt);
        camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
        delta_ints=zeros(1,PR_tries*2+1);
        
        %% Loops
        for j=1:big_it
            fprintf("big itteration %d\n",j);
            for i=1:PR_tries
                fprintf("PR try %d\n",i);
                %% Original solution
                camera_field_PR=GS_PR_max_mid(camera_int,N1);
                SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR))));
                SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                delta_ints(2*i-1)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                SLM_patterns(:,:,2*i-1)=SLM_patt_new;
                
                %% Conjugate solution
                camera_field_PR_conj=conj(camera_field_PR);
                SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR_conj))));
                SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                delta_ints(2*i)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                SLM_patterns(:,:,2*i)=SLM_patt_new;
            end
            [max_int,maxind]=max(delta_ints);
            SLM_patt=SLM_patterns(:,:,maxind);
            delta_ints(PR_tries*2+1)=max_int;
            SLM_patterns(:,:,PR_tries*2+1)=SLM_patt;
            camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
        end
        activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,1,forward_TM_show,N2,elems_c,elems_r);

    end
end