function [indices] = bad_indices(zscores,N)

k=1;

for i=1:N
    if zscores(i) <= 3
        indices(k)=i;
        k=k+1;
    end
end

if k==1
    indices=[];
end