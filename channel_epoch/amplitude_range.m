function zscore = amplitude_range(data,E,N)

for j=1:N
    amp_range(j) = max(data(j,:,E)) - min(data(j,:,E));
end

for j=1:N
    zscore(j) = abs(amp_range(j) - mean(amp_range))/std(amp_range);
end