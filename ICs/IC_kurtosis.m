function [zscore] = IC_kurtosis(ic_data,N)

%{
zscore= zeros(N,1);
num_time_frames=size(ic_data,2);
mu_four = zeros(N,1);
mu_two = zeros(N,1);
mu_two_square = zeros(N,1);
kurtosis_value = zeros(N,1);

for j=1:N
   mean_data = mean(ic_data(j,:));
   for i=1:num_time_frames
       mu_four(j) = mu_four(j)+ ( (ic_data(j,i)-mean_data).^4);
       mu_two(j) = mu_two(j)+ ( (ic_data(j,i)-mean_data).^2);
   end
   mu_four(j)= mu_four(j)/num_time_frames;
   mu_two(j)= mu_two(j)/num_time_frames;
   mu_two_square(j) = mu_two(j).^2;
   kurtosis_value(j) = (mu_four/mu_two_square) - 3;
end
%}

for i=1:N
    kurtosis_value(i)=kurt(ic_data(i,:));
end

for j=1:N
  zscore(j) = (abs(kurtosis_value(j) - median(kurtosis_value)))/MAD(kurtosis_value);
end