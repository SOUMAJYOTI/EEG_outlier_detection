function [zscore] = channels_variance(N,data)

zscore = zeros(N,1);
theta = zeros(N,1);

for j=1:N
   theta(j) = var(data(j,:));
end

median_theta = median(theta);
mean_theta = mean(theta);

MAD_theta = MAD(theta);
MEAN_theta = mean(abs(theta-mean(theta)));

for j=1:N
%  zscore(j) = (abs(theta(j) - median_theta))/MAD_theta;
  z_var(j,1) = (abs(theta(j) - median_theta))/MAD_theta;
end

for j=1:N
  zscore(j) = (abs(theta(j) - mean_theta))/MEAN_theta;
  z_var(j,2) = (abs(theta(j) - mean_theta))/MEAN_theta;
end


fprintf('\n');
%disp('THE ZSCORE VALUES FOR CHANNEL VARIANCES');
%disp('The first column displays MAD zscores and the second column displays the mean zscores for channel variances')
%display(z_var);


