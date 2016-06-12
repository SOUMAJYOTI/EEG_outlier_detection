function [rho]=spearman_coeff(x,y)
%Syntax: [rho]=spearman_coeff(x,y)
%__________________________
%
% Spearman's rank correalation coefficient.
%
% Written by Soumajyoti Sarkar

x_sort = sort(x);
count_map=zeros(size(unique(x_sort)));
curr_ptr=1;
count_curr=1;
k=1;
for i=2:size(x,2)
    if x_sort(i)==x_sort(curr_ptr)
        count_curr=count_curr+1;
    else
        count_map(k)=count_curr;
        k=k+1;
        curr_ptr=i;
        count_curr=1;
    end
end
count_map(k)=count_curr;
x_sort=unique(x_sort);
%disp(x_sort);
%disp(count_map);
k=1;
sum_pos=0.0;
i=1;
while i<=size(x,2)
    %disp(i)
    %disp(count_map(k))
    for j=0:count_map(k)-1
        sum_pos=sum_pos+j+i;
    end
    r=double(sum_pos)/double(count_map(k));
    for j=0:count_map(k)-1
        rank_x(i+j)=double(r);
    end
    i=i+count_map(k);
    k=k+1;
    sum_pos=0;
end
%rank_x(k)=r;
%disp(rank_x);

y_sort = sort(y);
count_map=zeros(size(unique(y_sort)));
curr_ptr=1;
count_curr=1;
k=1;
for i=2:size(y,2)
    if y_sort(i)==y_sort(curr_ptr)
        count_curr=count_curr+1;
    else
        count_map(k)=count_curr;
        k=k+1;
        curr_ptr=i;
        count_curr=1;
    end
end
count_map(k)=count_curr;
y_sort=unique(y_sort);
k=1;
sum_pos=0.0;
i=1;
while i<=size(y,2)
    %disp(i)
    %disp(count_map(k))
    for j=0:count_map(k)-1
        sum_pos=sum_pos+j+i;
    end
    r=double(sum_pos)/double(count_map(k));
    for j=0:count_map(k)-1
        rank_y(i+j)=double(r);
    end
    i=i+count_map(k);
    k=k+1;
    sum_pos=0;
end
%disp(rank_x)
%disp(rank_y)

%{
x_mean = mean(rank_x);
y_mean = mean(rank_y);
x_meanmat = repmat(x_mean, size(rank_x,2));
x_meanmat = x_meanmat(1,:);
y_meanmat = repmat(y_mean, size(rank_y,2));
y_meanmat = y_meanmat(1,:);

x_diff = rank_x-x_meanmat;
y_diff = rank_y-y_meanmat;
x_diff2 = power(x_diff,2);
y_diff2 = power(y_diff,2);
numer=0;

for i=1:size(rank_x,2)
    numer= numer+x_diff(i)*y_diff(i);
end
denom=sum(x_diff2)*sum(y_diff2);
denom=power(denom,0.5);
rho=numer/denom;
%}


end



