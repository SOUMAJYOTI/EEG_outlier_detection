function []= channels_statistics(data)

total_channels = size(data, 1);
total_frames = size(data, 2);

[z_score_channels_variance bad_channels_indices] = channels_variance(total_channels,data);
no_of_bad_channels=size(bad_channels_indices,2);
filtered_data = eeg_interp(data, bad_channels_indices, 'spherical');


end