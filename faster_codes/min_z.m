function [lengths] = min_z(list_properties)

y= ones(1,size(list_properties,2));
x= 3*ones(1,size(list_properties,2));

zs=list_properties-repmat(mean(list_properties,1),size(list_properties,1),1);
%display(zs);
zs=zs./repmat(std(zs,[],1),size(list_properties,1),1);
zs(isnan(zs))=0;
%all_l = abs(list_properties) > repmat(x,size(list_properties,1),1);
all_l = abs(zs) > repmat(x,size(list_properties,1),1);

lengths=zeros(size(all_l,1),1);

for i=1:size(all_l,1)
    for j=1:size(all_l,2)
        if all_l(i,j)==1
            lengths(i)=1;
        end
    end
end

%display(list_properties);
%display(zs);