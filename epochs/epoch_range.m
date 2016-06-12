function zscore = epoch_range(data,N,E)

for j=1:E

   for i=1:N

      range_epoch(j,i)=max(data(i,:,j))-min(data(i,:,j));
   end
   mean_range_epoch(j) = mean(range_epoch(j,:));
   
end

%display(mean_range_epoch);
for j=1:E
   zscore(j) = abs(mean_range_epoch(j) - median(mean_range_epoch))/MAD(mean_range_epoch);
end
%display(zscore);