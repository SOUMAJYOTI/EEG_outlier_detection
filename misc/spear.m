function [r]=spear(x,y)
%Syntax: [r,t,p]=spear(x,y)
%__________________________
%
% Spearman's rank correalation coefficient.
%
% r is the Spearman's rank correlation coefficient.
% t is the t-ratio of r.
% p is the corresponding p-value.
% x is the first data series (column).
% y is the second data series, a matrix which may contain one or multiple
%     columns.
%

% Find the data length
N = length(x);

% Get the ranks of x
%x=sort(x);
%y=sort(y);
R = crank(x)';
S = crank(y)';
   
%display(R);
%display(S);
% Calculate the correlation coefficient
r = 1-6*sum((R-S).^2)/N/(N^2-1);
    
end


function k=crank(x)
%Syntax: r=crank(x)
%__________________
%
% Assigns ranks on a data series x. 
%
% r is the vector of the ranks
% x is the data series. It must be sorted.
% 

u = unique(x);
[xs,z1] = sort(x);
[z1,z2] = sort(z1);
k = (1:length(x));
k=k(z2);

for i1=1:length(u)
    s=find(u(i1)==x); 
    k(s) = mean(k(s));
    
end

end




