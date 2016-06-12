function [eeg_channels]= channel_filter(EEG, rank_other)

EEGdata=EEG.data;
num_of_channels = size(EEGdata,1);
eeg_chans = 1:EEG.nbchan;

ref_chan=[];

measure = 1;
if strcmp(rank_other,'rank')==1
	% 1. rank Correlation
	zscore_channels_rankCorrelation = channels_rank_correlation_other_channels(num_of_channels, EEG.data);
	zs(:,measure) = zscore_channels_rankCorrelation;
else
	% 1. Pearson Correlation
	zscore_channels_Correlation = channels_correlation_other_channels(num_of_channels, EEG.data);
	zs(:,measure) = zscore_channels_Correlation;
end

measure = measure+1;
% 2. Variance
zscore_channels_variance = channels_variance(num_of_channels, EEG.data);
zs(:,measure) = zscore_channels_variance;

measure = measure+1;
%3. Hurst component
zscore_channels_hurstExponent = channels_Hurst_exponent(num_of_channels, EEG.data);
zs(:,measure) = zscore_channels_hurstExponent;

zs_channels = zs;
%display(zs_channels);
lengths = min_z_mod(zs);
chans_to_interp = eeg_chans(logical(lengths));
chans_to_interp = setdiff(chans_to_interp,ref_chan); % Ref chan may appear bad, but we shouldn't interpolate it!
display(chans_to_interp);

if(~isempty(chans_to_interp))
    channels_interpolated_eeg = eeg_interp(EEG, chans_to_interp);
end

eeg_channels=channels_interpolated_eeg;

end