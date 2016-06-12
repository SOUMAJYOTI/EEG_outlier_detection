function zscore = channels_rank_correlation_other_channels(N,data)

disp('Calculating Spearmans rank correlation coefficient...might take some time');

zscore = zeros(N,1);
theta = zeros(N,1);
%data=squeeze(data(1,:,:));
parfor j=1:N
   eta=zeros(N,1);
   for i=1:N
       %display(data(j,:));
      [corr] = spear(data(j,:),data(i,:));
%	disp(corr)
      if(isnan(corr)==1)
        eta(i) = 0;
      else
	%disp(corr);
        eta(i) = corr;
      end
   end
   theta(j) = mean (eta);
end

%display(theta)
median_theta = median(theta);
MAD_theta = MAD(theta);

%display(MAD_theta);
for j=1:N
   zscore(j) = (abs(theta(j) - median_theta))/MAD_theta;
%  zscore(j) = (abs(theta(j) - mean(theta)))/std(theta);
end