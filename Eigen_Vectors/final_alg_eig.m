%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations for finding eigen vectors of transmission matrix
clc;    close all;      clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init System..
num_elems=1;
[forward_TM_show,backward_TM_abs_orig,backward_TM_orig,N1,N2]=initialize_system_smooth();
O_support=get_circular_mask(N1,N1,N1/8);            %support for fourier field

element_location_cell{1} = create_element_matrix(N2,num_elems,[10 12 14 16 28 30 32 34]);
element_location_cell{2} = create_element_matrix(N2,num_elems,[46 47 48 49 50 51 52 53]);
element_location_cell{3} = create_element_matrix(N2,num_elems,[10 19 28 29 30]);
element_location_cell{4} = create_element_matrix(N2,num_elems,[10 18 73 81]);


for elems_i=1:4
    for kk=1:3
        element_location=element_location_cell{elems_i};
        elems=find(element_location>0);
        num_elems=numel(elems);
        num_iter=5;
        PR_tries=8;
        [elems_r,elems_c]=ind2sub([N2 N2],elems);
        forward_TM=forward_TM_show(elems,:);                     %
        backward_TM=backward_TM_orig(:,elems);                   %
        backward_TM_abs=backward_TM_abs_orig(:,elems);           %
        bt_tries=zeros(N1,N1,PR_tries*2+1);
        score_tries=zeros(1,PR_tries*2+1);
        %% First try
        bt_F=exp(1i*2*pi*rand(N1,N1)).*O_support;           %SLM pattern
        bt=ifftshift(ifft2(ifftshift(bt_F)));               %pattern on tissue
        bt=reshape(bt,N1^2,1);
        bt=bt/norm(bt);
        %% Loops
        fprintf('Number of camera shots: %d\n',num_iter*PR_tries*2);
        for i=1:num_iter
            fprintf('Itteration %d\n',i);
            at=abs((forward_TM*bt)).^2;                     % vector after tissue, only in fluor locations
            output_vector=backward_TM_abs*at;
            for j=1:PR_tries
                fprintf('PR try %d\n',j);
                bt = GS_PR(output_vector,N1);
                bt_conj=conj(bt);                               % valid solution to PR
                bt_tries(:,:,j*2-1)=bt;
                bt_tries(:,:,j*2)=bt_conj;
                
                %% Desicion between conj and orig
                bt_F=fftshift(fft2(fftshift(bt)));              %use bt on SLM
                bt_F=sign(bt_F).*O_support;                     %
                bt=ifftshift(ifft2(ifftshift(bt_F)));           %project bt
                bt=reshape(bt,N1^2,1);
                output_vec=backward_TM_abs*abs((forward_TM*bt)).^2;
                score_tries(j*2-1)=PR_score(output_vec);  %intensity in camera, using bt
                %         figure;imagesc(abs(reshape(output_vec,N1,N1)));title(j*2-1)
                
                bt_conj_F=fftshift(fft2(fftshift(bt_conj)));
                bt_conj_F=sign(bt_conj_F).*O_support;
                bt_conj=ifftshift(ifft2(ifftshift(bt_conj_F)));
                bt_conj=reshape(bt_conj,N1^2,1);
                output_vec_conj=backward_TM_abs*abs((forward_TM*bt_conj)).^2;
                score_tries(j*2)=PR_score(output_vec_conj);  %intensity in camera, using bt
                %         figure;imagesc(abs(reshape(output_vec_conj,N1,N1)));title(j*2)
                
            end
%             disp(score_tries)
            [max_score,maxind]=max(score_tries);
            bt=bt_tries(:,:,maxind);
            bt_tries(:,:,PR_tries*2+1)=bt;
            score_tries(PR_tries*2+1)=max_score;
            
            bt=reshape(bt,N1^2,1);
            a_show=forward_TM_show*(bt/norm(bt));
            a_show=reshape(a_show,N2,N2);
            if(i==num_iter)
                figure;imagesc(abs(a_show).^2);hold on;title([num2str(numel(elems)) ' elements in matrix']);colorbar
                plot(elems_c,elems_r,'*r');
            end
        end
        
        
    end
end