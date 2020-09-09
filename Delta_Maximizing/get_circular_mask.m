function [circ_mask]=get_circular_mask(M,N,mask_R,center_x,center_y)
if ~exist('center_x','var')
center_x=ceil(N/2);
end
if ~exist('center_y','var')
center_y=ceil(M/2);
end
[X,Y] = meshgrid(1:N,1:M);
distance = (X-center_x).^2+(Y-center_y).^2;
circ_mask = distance<mask_R^2;

