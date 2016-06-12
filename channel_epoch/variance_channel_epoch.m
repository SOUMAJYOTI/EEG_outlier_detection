function zscore = variance_channel_epoch(data,E,N)

%data=squeeze(data(1,:,:));
for j=1:N
    s(j) = var(data(j,:,E));
end


for j=1:N
    zscore(j) = abs(s(j) - mean(s))/std(s);   
end