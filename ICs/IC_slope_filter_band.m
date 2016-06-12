function zscore = IC_slope_filter_band(EEG,IC_data, f_LP1, f_LP2, N)

for u = 1:size(IC_data,1)
    [spectra(u,:) freqs] = pwelch(IC_data(u,:),[],[],(EEG.srate),EEG.srate);
end

for i=1:N
    mean_slope(i)= mean(diff(10*log10(spectra(i,find(freqs>=f_LP1,1):find(freqs<=f_LP2,1,'last')))));
end

for i=1:N
    zscore(i) = (abs(mean_slope(i) - median(mean_slope)))/MAD(mean_slope);
end

end