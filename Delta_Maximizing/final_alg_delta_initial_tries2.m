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
% for elems_i=1:numel(element_location_cell)
for elems_i=1:4
    for reps_i=1:3
        element_location=element_location_cell{elems_i};
        elems=find(element_location>0);
        num_elems=numel(elems);
        O_support=double(get_circular_mask(N1,N1,N1/8));                    %support for fourier field
        [elems_r,elems_c]=ind2sub([N2 N2],elems);                           %elemnt x,y
        forward_TM=forward_TM_show(elems,:);                                     %
        backward_TM=backward_TM_orig(:,elems);                                   %
        backward_TM_abs=backward_TM_abs_orig(:,elems);                           %
        PR_tries=2;
        big_it=3;
        delta_ints=zeros(1,PR_tries*2+1);
        SLM_patterns=zeros(N1,N1,PR_tries*2+1);
        alg_tries=20;
        final_delta_ints=zeros(1,alg_tries);
        final_SLM_patts=zeros(N1,N1,alg_tries);
        [X,Y]=meshgrid(1:N1);
        sigma=4;
        z=exp(-0.5/sigma.*((X-N1/2).^2+(Y-N1/2).^2));
        
        fprintf('Number of camera shots: %d\n',alg_tries*big_it*PR_tries*2);
        
        
        %% First try
        for alg_try_i=1:alg_tries
            fprintf("try new initialization %d out of %d\n",alg_try_i,alg_tries);
            SLM_patt=conv2(exp(1i*2*pi*rand(N1)),z,'same');  %first try - smooth rand phase
            SLM_patt=sign(SLM_patt).*O_support;
            camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
            delta_ints=zeros(1,PR_tries*2+1);
            
            %% Loops
            for j=1:big_it
                fprintf("big itteration %d our of %d\n",j,big_it);
                for i=1:PR_tries
                    fprintf("PR try %d out of %d\n",i,PR_tries);
                    %% Original solution
                    camera_field_PR=GS_PR_max_mid(camera_int,N1);
                    SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR))));
                    SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                    camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                    camera_int=normalize_var(camera_int);
                    delta_ints(2*i-1)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                    %                     delta_ints(2*i-1)=max(camera_int(:));
                    SLM_patterns(:,:,2*i-1)=SLM_patt_new;
                    %                     figure;imagesc(rot90(sqrt(camera_int),1));impixelinfo
                    
                    %% Conjugate solution
                    camera_field_PR_conj=conj(camera_field_PR);
                    SLM_patt_new=(fftshift(fft2(fftshift(camera_field_PR_conj))));
                    SLM_patt_new=sign(SLM_patt_new).*SLM_patt;
                    camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt_new,N1,0,forward_TM_show,N2,elems_c,elems_r);
                    camera_int=normalize_var(camera_int);
                    
                    delta_ints(2*i)=sum(sum(abs(camera_int(98:102,98:102)).^2));
                    %                     delta_ints(2*i)=max(camera_int(:));
                    SLM_patterns(:,:,2*i)=SLM_patt_new;
                    %                     figure;imagesc(rot90(sqrt(camera_int),1));impixelinfo
                    
                end
                [max_int,maxind]=max(delta_ints);
                SLM_patt=SLM_patterns(:,:,maxind);
                delta_ints(PR_tries*2+1)=max_int;
                SLM_patterns(:,:,PR_tries*2+1)=SLM_patt;
                camera_int = activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
                %                         figure;imagesc(rot90(sqrt(camera_int),1));impixelinfo
            end
            camera_int11=activate_vector_only_FFT(forward_TM,backward_TM,SLM_patt,N1,0,forward_TM_show,N2,elems_c,elems_r);
            final_SLM_patts(:,:,alg_try_i)=SLM_patt;
            final_delta_ints(alg_try_i)=max_int;
%             figure;imagesc(rot90(sqrt(camera_int11),3));impixelinfo
            
        end
        [~,max_delta_all_alg]=max(final_delta_ints);
        chosen_SLM_pat=final_SLM_patts(:,:,max_delta_all_alg);
        activate_vector_only_FFT(forward_TM,backward_TM,chosen_SLM_pat,N1,1,forward_TM_show,N2,elems_c,elems_r);
    end
end