function element_location = create_element_matrix(N,num_elems,locs)
element_location=zeros(N);
if(nargin==2 || isempty(locs))
    a=randperm(N^2,num_elems);
    element_location(a)=(1+rand(1,num_elems))/2;
    element_location=reshape(element_location,N^2,1);
elseif(nargin==3)
    element_location(locs)=(1+rand(1,numel(locs)))/2;
    element_location=reshape(element_location,N^2,1);
end