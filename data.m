function [noisy_data] = data(FRAMES, TRIALS, CHANNELS, noise_channel_percent,...
                             NOISE_CHAN_AMP)

% Inputs:
%  FRAMES - number of signal frames per each trial
%  TRIALS - number of simulated trials
%  SRATE - sampling rate of simulated signal
% Output:
%  synthetic_data - simulated EEG data; vector: 32 by frames*epochs containing concatenated trials


%% Load the RAW data
%filename = 'data_1.avr';     % This is the BESA generated data
%data_locs = 'data_1.elp';    % This is the BESA generated locations file

%% Convert the BESA data into EEGLAB format

%cnt = loadeep_avg(filename);    --------- %THIS PART IS NOT WORKING SINCE
                                           % IT REQURES EXTERNAL MEX FILES
                                           % THAT HAVE COMPILATION PROBLEMS
%display(cnt);

%%
%FRAMES=200;
%TRIALS=10;
SRATE=250;
% Generates synthetic data with the use of matlab script simulatedEEG.m
% Source - http://www.cs.bris.ac.uk/~rafal/phasereset/

%general parameters of the signal

%parameters of NE
NEAMP = -14;
NEFREQ = 5;
NEPOS = 150;
%parameters of PE
PEAMP = 11;
PEFREQ = 10;
PEPOS = 100;
%temporal jitter of peaks and amplitude of noise
TJITTER = 8;

disp ('Please wait while generating the data');

ne1 = NEAMP * peak (FRAMES, TRIALS, SRATE, NEFREQ, NEPOS, TJITTER);
pe1 = PEAMP * peak (FRAMES, TRIALS, SRATE, PEFREQ, PEPOS, TJITTER);
ne2 = 1.2*NEAMP * peak (FRAMES, TRIALS, SRATE, NEFREQ+2, NEPOS-100, TJITTER);

%display(size(ne1));
load dipole

dipole(32,:) = dipole(7,:);
if CHANNELS == 64
        dipole(33:64,:) = dipole(1:32,:);
end
if CHANNELS == 128
	dipole(33:64,:) = dipole(1:32,:);
        dipole(65:128,:) = dipole(1:64,:);
end

synthetic_data = dipole (:, 1) * ne1 + dipole (:, 4) * pe1 +dipole(:, 1) * ne2;
%disp(var(synthetic_data));
figure(1);plot(synthetic_data');
title('Original EEG data');

%save('noisy_datas/original_data', 'synthetic_data');
[FileName, PathName] = uiputfile('./noisy_datas/original_data/*.mat', 'Save original EEG data');
name = strcat(PathName,FileName);
save(name, 'synthetic_data');

disp('Generating noise for the channels');
fprintf('\n');

num_noisy_channels = round((noise_channel_percent*CHANNELS)/100);
perm_channels = randperm(CHANNELS);
noise_channels = perm_channels(1:num_noisy_channels);

%noise_channels = [2 12 15 20 22 25 26 30];
noise_channels = [sort(noise_channels)];
fprintf('\nThe %d percent noisy channels out of %d channels are:\n', noise_channel_percent, CHANNELS);
disp(noise_channels);
%noisy_data = synthetic_data;

for ch = 1:size(synthetic_data,1)
    %disp(sprintf('Generating noise for channel %d', ch));
    if(ismember(ch, noise_channels)==1)
	%disp(ch);
        noisy_data(ch,:) = synthetic_data(ch,:) + NOISE_CHAN_AMP * noise (FRAMES, TRIALS, SRATE);
    else
        noisy_data(ch,:) = synthetic_data(ch,:);
    end
end
%disp(synthetic_data(6,:));
figure(2), plot(noisy_data');
title('EEG data after adding noise');
%name= strcat('noisy_datas/noisy_channel_data_',num2str(noise_channel_percent),'_percent_channel', num2str(NOISE_CHAN_AMP),'_amplitude');
[FileName, PathName] = uiputfile('./noisy_datas/noisy_channel_data/*.mat', 'Save noisy channel EEG data');
name = strcat(PathName,FileName);
save(name, 'noisy_data');

%{
save('data_model/noisy_epoch_data', 'noisy_epoch_data');

%}
%%
%figure(2);plot(noisy_data');
% fileID = fopen('ascii_eeg.txt','w');
% [L,W]=size(synthetic_data);
% for i=1:L
%     for j=1:W
%         fprintf(fileID, '%f ',synthetic_data(i,j));
%     end
%     fprintf(fileID,'\n');
% end

%%

