function zscore = IC_median_gradient(IC_data,N)

slope_spatial_component = zeros(N,1);
%IC_time_derivative = zeros(IC_time,1);
zscore = zeros(N,1);

for j=1:N
    slope_spatial_component(j)= median(diff(IC_data(j,:)));
end

for i=1:N
    zscore(i) = (abs(slope_spatial_component(i) - median(slope_spatial_component)))/MAD(slope_spatial_component);
end

end