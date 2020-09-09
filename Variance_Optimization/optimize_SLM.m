function best_SLM_patt = optimize_SLM(forward_TM,N1,backward_TM_abs,O_support,forward_TM_show,N2)

modif_eps=0.825;
modif_eps_power=0.5;

taking_best=5;
random_tries=5;
first_try_attempts=taking_best*(random_tries+1);
modif_loops=200;

h_variances=zeros(1,first_try_attempts);
taken_vectors=zeros(N1,N1,taking_best);
input_vector_SLMs=zeros(N1,N1,first_try_attempts);
a_plots=zeros(size(forward_TM,1),modif_loops);
scores=zeros(1,modif_loops);
fprintf('number of camera shots: %d\n',first_try_attempts*modif_loops*taking_best*random_tries);


%% First Guesses
for ind=1:first_try_attempts
 
    b_in_SLM=exp(1i*2*pi*rand(N1,N1)).*O_support;
    b_in=fftshift(fft2(fftshift(b_in_SLM)));
    input_vector_SLMs(:,:,ind)=b_in_SLM;
    [b_out] = activate_vector(forward_TM,backward_TM_abs,b_in,0,forward_TM_show,N2);
    [h_variances(ind)]=vec_score(b_out);
end

[h_vars,prev_indeces]=sort(h_variances,'descend','MissingPlacement','last');
taken_vectors(:,:,1:taking_best)=input_vector_SLMs(:,:,prev_indeces(1:taking_best));
h_variances(1:(1+random_tries):taking_best*(random_tries+1))=h_vars(1:taking_best);
input_vector_SLMs(:,:,1:(1+random_tries):taking_best*(random_tries+1))=taken_vectors(:,:,1:taking_best);

%% converging loops
for loop=1:modif_loops
    if(mod(loop,100)==0)
        fprintf('loop: %d out of %d\n',loop,modif_loops);
    end
    SLM_diff_arg=(modif_eps*(loop^-modif_eps_power));
    ind=1;
    for good_input_vec_i=1:taking_best
        %% debug
        if(good_input_vec_i==1)
            b_in_SLM=taken_vectors(:,:,1);
            b_in=fftshift(fft2(fftshift(b_in_SLM)));
            [b_out,at] = activate_vector(forward_TM,backward_TM_abs,b_in,1,forward_TM_show,N2);
            a_plots(:,loop)=at;
            scores(loop)=vec_score(b_out);
            figure(123);
            subplot(2,1,1);plot(a_plots')
            subplot(2,1,2);plot(scores)
        end
        
        %% continue
        b_in_SLM=taken_vectors(:,:,good_input_vec_i);
        ind=ind+1;
        for small_mod_i=1:random_tries
            SLM_diff=(exp(1i*2*pi*SLM_diff_arg*(rand(N1,N1)-0.5)));
            moded_b_slm=b_in_SLM.*SLM_diff;
            input_vector_SLMs(:,:,ind)=moded_b_slm;
            b_in=fftshift(fft2(fftshift(moded_b_slm)));
            [b_out] = activate_vector(forward_TM,backward_TM_abs,b_in,0,forward_TM_show,N2);
            [h_variances(ind)]=vec_score(b_out);
            ind=ind+1;
        end
    end
    [h_vars,prev_indeces]=sort(h_variances,'descend','MissingPlacement','last');
    taken_vectors(:,:,1:taking_best)=input_vector_SLMs(:,:,prev_indeces(1:taking_best));
    h_variances(1:(1+random_tries):taking_best*(random_tries+1))=h_vars(1:taking_best);
    input_vector_SLMs(:,:,1:(1+random_tries):taking_best*(random_tries+1))=taken_vectors(:,:,1:taking_best);
    
    
end
best_SLM_patt=taken_vectors(:,:,1);
end