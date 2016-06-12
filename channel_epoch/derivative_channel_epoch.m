function zscore = derivative_channel_epoch(data,E,N)

%data=squeeze(data(1,:,:));
for j=1:N
  d(j) = median(diff(data(j,:,E)));
end

for j=1:N
    zscore(j) = abs(d(j) - median(d))/MAD(d);
end
