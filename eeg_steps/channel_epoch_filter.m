function [chans_to_interp]= channel_epoch_filter(EEG, epoch)

eeg_chans=1:size(EEG.data,1);
ref_chan=[];
epoch_data = EEG.data;
num_of_channels = size(epoch_data,1);
num_of_time_points = size(epoch_data, 2);
num_of_epochs = size(epoch_data,3);

measure=1;
zscore_epoch_variance = variance_channel_epoch(epoch_data, epoch ,num_of_channels);
zs(:,measure) = zscore_epoch_variance;
measure=measure+1;

%zscore_epoch_mean = mean_channel_epoch(epoch_data, epoch, num_of_channels);
%zs(:,2) = zscore_epoch_mean;

zscore_epoch_derivative = derivative_channel_epoch(epoch_data, epoch, num_of_channels);
zs(:,measure) = zscore_epoch_derivative;
measure=measure+1;

zscore_epoch_amplitude = amplitude_range(epoch_data, epoch, num_of_channels);
zs(:,measure) = zscore_epoch_amplitude;
measure=measure+1;

%disp(zs);
lengths = min_z_mod(zs);
chans_to_interp = eeg_chans(logical(lengths));

%{
chans_to_interp = setdiff(chans_to_interp,ref_chan); % Ref chan may appear bad, but we shouldn't interpolate it!
fprintf('The channels to interpolate are:\n');
disp(chans_to_interp);

if(~isempty(chans_to_interp))
    EEG = eeg_interp(EEG, chans_to_interp);
end

channel_epoch_structure = EEG;
%eeg_channels_data = channels_interpolated;
%}

end