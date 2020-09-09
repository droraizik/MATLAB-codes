function [assumed_field_final] = GS_PR(intens_image,N)
    %% Algorithm
    num_iter=500;
    intens_image=sqrt(reshape(intens_image,N,N));
    O_support=get_circular_mask(N,N,N/8);
    first_try=(intens_image).*exp(1i*2*pi*rand(N,N));
    assumed_F=fftshift(fft2(fftshift(first_try))).*O_support;
    for i=1:num_iter
        assumed_field=ifftshift(ifft2(ifftshift(assumed_F)));
        assumed_field=(intens_image).*sign(assumed_field);
        assumed_F=fftshift(fft2(fftshift(assumed_field)));
        assumed_F=assumed_F.*O_support;
    end
    assumed_field_final=ifftshift(ifft2(ifftshift(assumed_F)));
end
