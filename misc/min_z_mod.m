function [lengths] = min_z_mod(zs)

x= 3*ones(1,size(zs,2));

zs(isnan(zs))=0;
all_l = abs(zs) > repmat(x,size(zs,1),1);
lengths=zeros(size(all_l,1),1);

for i=1:size(all_l,1)
    for j=1:size(all_l,2)
        if all_l(i,j)==1
            lengths(i)=1;
        end
    end
end
%display(lengths);