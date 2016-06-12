function zscore = channels_correlation_other_channels(N,data)

zscore=zeros(N,1);
theta = zeros(N,1);

for j=1:N
   eta=zeros(N,1);
   for i=1:N
      cor = abs(corrcoef(data(j,:),data(i,:)));
      if(isnan(cor(1,2))==1)
        eta(i) = 0;
      else
        eta(i) = cor(1,2);
      end
   end
   theta(j) = mean (eta);
end

median_theta = median(theta);
MAD_theta = MAD(theta);

%display(theta);
for j=1:N
  zscore(j) = (abs(theta(j) - median_theta))/MAD_theta;
 % zscore(j) = (abs(theta(j) - mean(theta)))/std(theta);
end
%display(zscore);