function zscore = epoch_mean_median(data, E)

eeg_chans=1:32;
medians = median(data(eeg_chans,:),2);

%  Epoch's mean deviation from channel means.
for u = 1:size(data,3)
	mean_median_epoch(u) = mean(abs(squeeze(mean(data(eeg_chans,:,u),2)) - medians));
end

%{
for j=1:E
   for i=1:N
      mean_median_epoch(j,i)=mean(data(i,:,j)) ;
   end
   mean_mean_median_epoch(j) = mean(mean_median_epoch(j,:));
   %display(mean_mean_median_epoch(j));
end
%}

for j=1:E
   zscore(j) = abs(mean_median_epoch(j) - median(mean_median_epoch))/MAD(mean_median_epoch);
end
