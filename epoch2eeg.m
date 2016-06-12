function [concat_data]=epoch2eeg(data)

frames = size(data,2);
k=0;
for i=1:size(data,1)
    k=0;
    for j=1:size(data,3)
        concat_data(i,k+1:k+frames)=data(i,:,j);
        k=k+frames;
    end
end
end