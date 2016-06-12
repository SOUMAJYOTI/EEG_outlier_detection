function [filtered_epoch epoch_noisy]= epoch_filter(EEG)

epoch_data = EEG.data;
num_of_channels = size(epoch_data,1);
num_of_epochs = size(epoch_data,3);
measure = 1;

zscore_epoch_mean_median = epoch_mean_median(epoch_data, num_of_epochs);
zs(:,measure) = zscore_epoch_mean_median;
measure = measure+1;

zscore_epoch_variance = epoch_variance(epoch_data,num_of_channels,num_of_epochs);
zs(:,measure) = zscore_epoch_variance;
measure=measure+1;

zscore_epoch_range = epoch_range(epoch_data,num_of_channels,num_of_epochs);
zs(:,measure) = zscore_epoch_range;
measure=measure+1;

display(zs);
lengths = min_z_mod(zs);
%disp(lengths);
epoch_noisy=[];
k=1;
for i=1:num_of_epochs
    if lengths(i)==1
        epoch_noisy(k)=i;
        k=k+1;
    end
end

fprintf('The noisy epochs detected by our method are\n\n');
if ~isempty(epoch_noisy)
    display(epoch_noisy-1);
else
    fprintf('\tNo noisy epochs detected!!\n');
end
    
k=1;
for i=1:size(epoch_data,3)
    if ismember(i,epoch_noisy)==0
        filtered_epoch(:,:,k) = epoch_data(:,:,i); 
        k=k+1;
    end
end

%EEG.data = filtered_epoch;
%EEG=pop_rejepoch(EEG, find(lengths),0);


%epoch_filtered_structure = EEG;
%eeg_channels_data = channels_interpolated;

end