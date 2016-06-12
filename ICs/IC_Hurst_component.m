function zscore = IC_Hurst_component(IC_data, N)

for j=1:N
   theta(j) = estimate_hurst_exponent(IC_data(j,:));
end

for i=1:N
    zscore(i) = (abs(theta(i) - median(theta)))/MAD(theta);
end

end