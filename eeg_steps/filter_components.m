function [ica] = filter_components(icaEEG, EEG, ica_data, lopass_freq)

num_components_pca=size(ica_data,1);

blink_chans = [4 11 22 29];
for i=1:4
    blink_channel_data(i,:)= EEG.data(blink_chans(i),:);
end

zscore_ica_correlation_EOG = IC_correlation_EOG_channels(ica_data, blink_channel_data);
zs(:,1) = zscore_ica_correlation_EOG;
%display(zs);

zscore_ica_kurtosis = IC_kurtosis(ica_data, num_components_pca);
zs(:,2) = zscore_ica_kurtosis;


zscore_ica_median_gradient = IC_median_gradient(ica_data, num_components_pca);
zs(:,3) = zscore_ica_median_gradient;

f_LP1=lopass_freq-5;
f_LP2=lopass_freq+5;
zscore_ica_slope_filter_band = IC_slope_filter_band(EEG, ica_data, f_LP1, f_LP2, num_components_pca);
zs(:,4) = zscore_ica_slope_filter_band;

zscore_ica_Hurst = IC_Hurst_component(ica_data, num_components_pca);
%display(zscore_ica_Hurst);
zs(:,5) = zscore_ica_Hurst;

[lengths] = min_z_mod(zs);

%display(zs);
% Reject
if ~isempty(find(lengths,1))
    EEGdata = pop_subcomp_mod(icaEEG, find(lengths), 0);
end

ica =EEGdata;

end