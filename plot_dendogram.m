function plot_dendogram(epoched_data)

num_epochs = size(epoched_data,3);
for i=1:size(epoched_data,3)
    for j=1:size(epoched_data,3)
        for k=1:size(epoched_data,3)
            if k~=i && k~=j
                dist_matrix(i,j) = abs((norm(epoched_data(:,:,i)-epoched_data(:,:,k)) - norm(epoched_data(:,:,j)-epoched_data(:,:,k))))/(num_epochs-2);
            end
        end
        
    end
    
end

disp((dist_matrix));
Z=linkage(dist_matrix);
disp(Z);
figure();
dendrogram(Z);
end