function zscore = epoch_variance(data,N,E)

%data=squeeze(data(1,:,:));
%display(data(1,:,1));
for j=1:E
   for i=1:N
      variance_epoch(j,i)=var(data(i,:,j));
   end
   mean_variance_epoch(j) = mean(variance_epoch(j,:));
   %display(mean_variance_epoch);
end

for j=1:E
   zscore(j) = abs(mean_variance_epoch(j) - median(mean_variance_epoch))/MAD(mean_variance_epoch);
end