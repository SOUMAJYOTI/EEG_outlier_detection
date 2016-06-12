function zscore = range_channel_epoch(data,E,N)

%data=squeeze(data(1,:,:));
for j=1:E
  for i=1:N
  r(j,i) = max(data(i,(15*(j-1)+1):15*j))-min(data(i,(15*(j-1)+1):15*j));
  end
end

for j=1:E
  for i=1:N
    zscore(j) = abs(r(j,i) - median(r(j,:)))/MAD(r(j,:));
  end
end