function zscore = mean_channel_epoch(data,E,N)

%data=squeeze(data(1,:,:));
for j=1:N
      m(j)=mean(data(j,:,E))- median(median(data(j,:,:)));
end

for j=1:N
    zscore(j) = abs(m(j) - mean(m))/std(m);
end
%list_properties(:,measure)=abs(mean(EEG.data(eeg_chans,:,epoch_num),2)-mean(EEG.data(eeg_chans,:),2));