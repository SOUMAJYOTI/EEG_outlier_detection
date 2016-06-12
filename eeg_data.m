function eeg_data(FRAMES, TRIALS, CHANNELS, NOISE_CHANNEL_PERCENT, NOISE_CH_AMP,...
                    NOISE_EPOCH_PERCENT, NOISE_EP_AMP)

                
clc;

SRATE=250; % this is fixed

%parameters of NE
NEAMP = -2; %-3
NEFREQ = 2;
NEPOS = 150;
%parameters of PE
PEAMP = 2; %3
PEFREQ = 2;
PEPOS = 60;
%temporal jitter of peaks and amplitude of noise
TJITTER = 8;

disp ('Please wait while generating the data');
load dipole

dipole(32,:) = dipole(7,:);
if CHANNELS == 64
        dipole(33:64,:) = dipole(1:32,:);
end
if CHANNELS == 128
	dipole(33:64,:) = dipole(1:32,:);
    dipole(65:128,:) = dipole(1:64,:);
end

TRIALS=TRIALS+1;  % the 1 added is for the baseline trial at the beginning
             
ne4 = NEAMP * peak (FRAMES, 1, SRATE, NEFREQ+2, 100, TJITTER);
ne5 = 0.7*PEAMP * peak (FRAMES, 1, SRATE, PEFREQ+2, 150, TJITTER);
epoch_noise = dipole (:, 1) * ne4 +dipole(:, 1) * ne5;
%epochN = epoch_noise(:,FRAMES+1:(FRAMES*2));  %epochN is used as the epoch noisy data to be added to the EEG


%%%%%%%%%%%%%%%% RANDOMLY SELECT NOISY CHANNELS %%%%%%%%%%%%%%%%%%%%%%%%
num_noisy_channels = round((NOISE_CHANNEL_PERCENT*CHANNELS)/100);
perm_channels = randperm(CHANNELS);
noise_channels = perm_channels(1:num_noisy_channels);
fprintf('\nThe %d percent noisy channels out of %d channels are:\n', NOISE_CHANNEL_PERCENT, CHANNELS);
disp(noise_channels);

%%%%%%%%%%%%%%%% RANDOMLY SELECT NOISY EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%
num_noisy_epochs = round((NOISE_EPOCH_PERCENT*TRIALS)/100);
perm_epochs = randperm(TRIALS+1);
noise_epochs = perm_epochs(1:num_noisy_epochs);
fprintf('\nThe %d percent noisy epochs out of %d epochs are:\n', NOISE_EPOCH_PERCENT, TRIALS);
disp(noise_epochs-1)

num_noises= mod(round(rand(1)*50), 5)+1; %randomizes the number of trials/places to introduce undulations in original EEG
%disp(num_noises);
num_noises_channels= mod(round(rand(1)*50), TRIALS-1)+1; %randomizes the number of trials/places to introduce undulations in original EEG
%disp(num_noises_channels);

count_5=0;
for i=1:num_noises
    rand_trial = mod(round(rand(1)*50),TRIALS-1)+1;
    
    % place stores the locations of the undulations in original EEG
    place(i,1)= FRAMES*rand_trial+15;
    place(i,2)= FRAMES*rand_trial+70;
    place(i,3)= FRAMES*rand_trial+135;
%    disp(rand_trial);
end

for i=1:num_noises_channels
    rand_trial = mod(round(rand(1)*50),TRIALS-1)+1;
    
    % place stores the locations of the channel noises
    place_channels(i,1)= FRAMES*rand_trial+150;
    place_channels(i,2)= FRAMES*rand_trial+60;
    place_channels(i,3)= FRAMES*rand_trial+100;
    place_channels(i,4)= FRAMES*rand_trial+180;
    %disp(rand_trial);
end

for i=1:CHANNELS
    mynoise = noise (FRAMES*TRIALS, 1, SRATE);
    mysignal(i,:)=zeros(1,FRAMES*TRIALS);
    for j=1:num_noises
        mypeak1 = peak (FRAMES*TRIALS, 1, SRATE, 3, place(j,1));
        mypeak2 = peak (FRAMES*TRIALS, 1, SRATE, 5, place(j,2));
        mypeak3 = peak (FRAMES*TRIALS, 1, SRATE, 0.8, place(j,3));
    
        a=round(rand(1)*10);
        if mod(a,7)==0
            sign=0;
        elseif mod(a,2)==1
            sign=-1;
        else
            sign=1;
        end
        mysignal(i,:) = mysignal(i,:)+ sign * mypeak1 + 0.6*sign*mypeak2 + sign*mypeak3;
    end
    mysignal(i,:) = mysignal(i,:)+4* mynoise;
    
    synthetic_data=zeros(CHANNELS,FRAMES*TRIALS);
    for j=1:num_noises_channels
        ne1 = NEAMP * peak (FRAMES*TRIALS, 1, SRATE, NEFREQ, place_channels(j,1), TJITTER);
        pe1 = PEAMP * peak (FRAMES*TRIALS, 1, SRATE, PEFREQ, place_channels(j,2), TJITTER);
        ne2 = 0.7*NEAMP * peak (FRAMES*TRIALS, 1, SRATE, NEFREQ-1, place_channels(j,3), TJITTER);
        ne3 = 0.7*NEAMP * peak (FRAMES*TRIALS, 1, SRATE, NEFREQ-1, place_channels(j,4), TJITTER);

        synthetic_data = synthetic_data + dipole (:, 1) * ne1 + dipole (:, 4) * pe1 +...
                 dipole(:, 1) * ne2 + dipole(:, 1) * ne3 ;      % the synthetic data is used as the channel noise, added to the original EEG
    end
%    plot(synthetic_data');
    a=round(rand(1)*10);
        if mod(a,7)==0
            sign=0;
        elseif mod(a,2)==1
            sign=-1;
        else
            sign=1;
        end
    if(ismember(i, noise_channels)==1)
        noisy_signal(i,:) =   mysignal(i,:) + sign*NOISE_CH_AMP*synthetic_data(i,:);
    else
        noisy_signal(i,:) = mysignal(i,:);
    end
   
    
end

x = [1:(TRIALS)*FRAMES];
x = x-FRAMES;
figure(1), plot(x,mysignal'), axis([-FRAMES ((TRIALS-1)*FRAMES) -15 15]);
title('Original EEG data');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

save_data(1,:,:) = mysignal;
%{
[FileName, PathName] = uiputfile('./noisy_datas/*.mat', 'Save original EEG data');
name = strcat(PathName,FileName);
save(name, 'mysignal');
%}

x = [1:(TRIALS)*FRAMES];
x = x-FRAMES;
figure(2), plot(x,noisy_signal'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
title('EEG after introducing the dipoles in the channels');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
%{
[FileName, PathName] = uiputfile('./noisy_datas/*.mat', 'Save noisy EEG data');
name = strcat(PathName,FileName);
save(name, 'noisy_signal');
%}

%figure(4),plot( NOISE_EP_AMP*epochN')

for tr=1:TRIALS
        if(ismember(tr, noise_epochs )==1)
            noisy_signal(:,(tr-1)*FRAMES+1:(tr)*FRAMES) =   noisy_signal(:,(tr-1)*FRAMES+1:(tr)*FRAMES) + NOISE_EP_AMP*epoch_noise; 
        end
end

save_data(2,:,:) = noisy_signal(:,:);
figure(3), plot(x,noisy_signal'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
title('EEG after introducing the dipoles in the epochs');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
[FileName, PathName] = uiputfile('./noisy_datas/*.mat', 'Save EEG datas');
name = strcat(PathName,FileName);
save(name, 'save_data');
%}
%}

end