function [zscore] = IC_correlation_EOG_channels(ic_data,eog_data)

N=size(ic_data,1);
zscore= zeros(N,1);
%num_time_frames=size(ic_data,2);
correlation_EOG_values = zeros(N,4);
correlation_EOG_max_value = zeros(N,1);

%plot(eog_data(2,:));
for j=1:size(ic_data,1)
    for i=1:4
     corr=corrcoef(ic_data(j,:),eog_data(i,:));
     correlation_EOG_values(j,i)=abs(corr(1,2));
    end
    correlation_EOG_max_value(j)=max(correlation_EOG_values(j,:));
end

%display(correlation_EOG_max_value);
for j=1:size(ic_data,1)
  zscore(j) = (abs(correlation_EOG_max_value(j) - median(correlation_EOG_max_value)))/MAD(correlation_EOG_max_value);
end