function zscore = channels_Hurst_exponent(N,data)

zscore= zeros(N,1);
theta = zeros(N,1);

for j=1:N
   theta(j) = estimate_hurst_exponent(data(j,:));
   if(isnan(theta(j))==1)
       theta(j)=0;
   end
end

MAD_theta = median(abs(theta-0.7));

for j=1:N
  zscore(j) = (abs(theta(j) - 0.7))/MAD_theta;
end