function main_script(TRIALS, CHANNELS, NAME_DATA, isRank)

% main_script() - Processes all the stages in the EEG artifact rejection
%
% Usage:
%  >> main_script(5, 64, 'noisy_datas/1.mat'); 
%
% Inputs:
%   TRIALS- number of epochs
%   CHANNELS - number of channels

%clear all;

clc
SRATE =250;
fig_count=1;
FRAMES = 200;

TRIALS = TRIALS+1;
load(NAME_DATA);

if nargin<4 
    isRank=1;
end

original_signal(:,:) = save_data(1,:,:);
noisy_signal(:,:) = save_data(2,:,:);

%disp(size(original_signal));
x = [1:(TRIALS)*FRAMES];
x = x-FRAMES;
figure(fig_count), plot(x,original_signal'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
fig_count=fig_count+1;
title('Raw EEG');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

x = [1:(TRIALS)*FRAMES];
x = x-FRAMES;
figure(fig_count), plot(x,noisy_signal'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
fig_count=fig_count+1;
title('Original noisy EEG');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

rawInputEEG = pop_importdata('data',noisy_signal);
        
if CHANNELS==32
    electrode_filename='data_sample/32.ced';
elseif CHANNELS==64
    electrode_filename='data_sample/64.ced';
else
    electrode_filename='data_sample/128.ced';
end

rawInputEEG.chanlocs=readlocs(electrode_filename);

%disp(size(noisy_signal));
%------------------------------------------------------------

    % Write EEG ascii data to file
    %fileID = fopen('d:\work\code\data_sample\ascii_eeg.txt','w');
    %{
    [L,W]=size(rawInputEEGdata);
    for i=1:L
        for j=1:W
            fprintf(fileID, '%f ',rawInputEEGadata(i,j));
        end
        fprintf(fileID,'\n');
    end
    %}
%------------------------------------------------------------

    %Plot the EEG data
    %eegplot(rawInputEEGdata,'eloc_file',electrode_filename);
%------------------------------------------------------------
EEG=rawInputEEG;
%Reference data to average
EEG=pop_reref(EEG, [], 'keepref', 'on');
epochEEG=EEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpf_on=0;
lpf_on=0;
notch_on=0;

if hpf_on ==0
    fprintf('\nHigh Pass Filtering off\n');
end
    
if lpf_on == 0
    fprintf('\nLow Pass Filtering off\n');
end
    
if notch_on ==0
    fprintf('\nNotch Pass Filtering off\n');
end
    
do_hipass=hpf_on;
do_lopass=lpf_on;
do_notch=notch_on;

hpf_freq=1;
hpf_bandwidth=0.50;
hpf_ripple=0.05;
hpf_attenuation=100;

lpf_freq=95;
lpf_bandwidth=2.50;
lpf_ripple=0.01;
lpf_attenuation=40;

notch_freq=50;
notch_bandwidth1=3;
notch_bandwidth2=1;
notch_ripple=0.05;
notch_attenuation=80;
   
% High Pass filtering
%----------------------
if do_hipass
    w_h=hpf_freq;
    t_h=hpf_bandwidth;
    r_h=hpf_ripple;
    a_h=hpf_attenuation;
    [m, wtpass, wtstop] = pop_firpmord([w_h-(t_h) w_h+(t_h)], [0 1], [10^(-1*abs(a_h)/20) (10^(r_h/20)-1)/(10^(r_h/20)+1)], EEG.srate);
    if mod(m,2);m=m+1;end;
    EEG = pop_firpm(EEG, 'fcutoff', w_h, 'ftrans', t_h, 'ftype', 'highpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
end

% Low Pass filtering
%----------------------
if do_lopass
    w_l=lpf_freq;
    t_l=lpf_bandwidth;
    r_l=lpf_ripple;
    a_l=lpf_attenuation;
  
    [m, wtpass, wtstop] = pop_firpmord([w_l-(t_l) w_l+(t_l)], [1 0], [(10^(r_l/20)-1)/(10^(r_l/20)+1) 10^(-1*abs(a_l)/20)], EEG.srate);
    if mod(m,2);m=m+1;end;
    EEG = pop_firpm(EEG, 'fcutoff', w_l, 'ftrans', t_l, 'ftype', 'lowpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
end

% NOTCH FILTER
%----------------------
if do_notch
    for n=1:length(o.filter_options.notch_freq)
        w_n=[notch_freq(n)-notch_bandwidth1/2 notch_freq(n)+notch_bandwidth1/2];
        t_n=notch_bandwidth2;
        r_n=notch_ripple;
        a_n=notch_attenuation;
        [m, wtpass, wtstop] = pop_firpmord([w_n(1)-(t_n) w_n(1)+(t_n) w_n(2)-(t_n) w_n(2)+(t_n)], [0 1 0], [10^(-1*abs(a_n)/20) (10^(r_n/20)-1)/(10^(r_n/20)+1) 10^(-1*abs(a_n)/20)], EEG.srate);
        if mod(m,2);m=m+1;end;
        EEG = pop_firpm(EEG, 'fcutoff', w_n, 'ftrans', t_n, 'ftype', 'bandstop', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
    end
end


tempEEG = EEG;
%-------------------------------------------------------------
%OUR MODIFIED FASTER METHOD FOR CHANNEL FILTERING
%-------------------------------------------------------------

fprintf('\n');
disp('STAGE1: OUR METHOD FOR CHANNEL FILTERING');
%Each Channel examination
%figure,plot(refEEG.data(7,:));
%eegplot(refEEG.data,'eloc_file',electrode_filename);

if isRank
    fprintf('\nWorking with Rank correlation\n' );
    [eeg_channels]= channel_filter(EEG,'rank');
    EEG = eeg_channels;
    title('Stage 1:Channel filtered EEG data using Rank correlation');
    x=1:size(EEG.data,2);
    x=x-FRAMES;
    figure(fig_count),plot(x,EEG.data'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
    fig_count=fig_count+1;
    %name = strcat('saved_data_model/chanFiltered_RankCorrelation_', num2str(noise_channel_perc),'_percent_noise_',num2str(noise_amp),'_noise_amp');
    %save(name, 'EEG');
    [FileName, PathName] = uiputfile('./saved_data_outputs/*.mat', 'Save channel filtered EEG data using rank correlation');
    name = strcat(PathName,FileName);
    EEGdata = EEG.data;
    save(name, 'EEGdata');
    epochEEG = EEG;
end

%figure(fig_count),plot(EEG.data');
%title('Channel filtered EEG data using Rank correlation');
fprintf('\nWorking with Pearson correlation \n');
[eeg_channels]= channel_filter(tempEEG,'other');
EEG = eeg_channels;
%name = strcat('saved_data_model/chanFiltered_PearsonCorrelation_',num2str(noise_channel_perc),'_percent_noise_',num2str(noise_amp),'_noise_amp');
%save(name, 'EEG');
x=1:size(EEG.data,2);
x=x-FRAMES;
figure(fig_count),plot(x,EEG.data'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
title('Stage 1:Channel filtered EEG using Pearson correlation');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;
[FileName, PathName] = uiputfile('./saved_data_outputs/*.mat', 'Save channel filtered EEG data using Pearson correlation');
name = strcat(PathName,FileName);
EEGdata = EEG.data;
save(name, 'EEGdata');

%----------------------------------------------------
%   BASELINE CORRECTION BY DEDUCTING THE AVERAGE
%----------------------------------------------------

data_average_baseline = mean(EEG.data(:,1:FRAMES));
data_correct = repmat(data_average_baseline, CHANNELS, 1);

%disp(size(data_correct));
for i=1:TRIALS
    EEG.data(:, (i-1)*FRAMES+1:i*FRAMES) = EEG.data(:, (i-1)*FRAMES+1:i*FRAMES) - data_correct;
end
x=1:size(EEG.data,2);
x=x-FRAMES;
figure(fig_count),plot(data_average_baseline')%, axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
title('Baseline average after Stage 1');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;

%--------------------------------------------------------------
% EPOCH FILTERING METHOD USING OUR MODIFIED FASTER
%--------------------------------------------------------------

fprintf('\n\n');
disp('STAGE 2: Starting Epoch filtering stage using our modified FASTER method');
fprintf('\n');

epdata = EEG.data;
epoched_data = epoch_data(FRAMES, TRIALS, epdata);

epochEEG.data = epoched_data;
epochEEG.pnts = size(epoched_data,2);
epochEEG.srate = SRATE;
epochEEG.trials = size(epoched_data,3);

tempEpoch= epochEEG;  % Saving a copy of channel filtered EEG for epoch filtering using FASTER
%disp(epochEEG);

%disp('The zscores from epoch properties calculated by our method are:');
fprintf('\n');
%fprintf('\tMEAN_MEDIAN  VARIANCE  RANGE\n');
[filtered_epoch epoch_noisy] = epoch_filter(epochEEG);

filterEpoch = epoch2eeg(filtered_epoch);

epochEEG.data = filtered_epoch;
epochEEG.pnts = size(filtered_epoch,2);
epochEEG.srate = SRATE;
epochEEG.trials = size(filtered_epoch,3);
TRIALS = size(filtered_epoch,3);

x=1:size(filterEpoch,2);
x=x-FRAMES;
figure(fig_count),plot(x,filterEpoch'), axis([-FRAMES size(filterEpoch,2)-FRAMES -15 15]);
fig_count=fig_count+1;
title('Stage2:Epoch filtered EEG data using our method');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

[FileName, PathName] = uiputfile('./saved_data_outputs/*.mat', 'Save epoch filtered EEG data using our method');
name = strcat(PathName,FileName);
save(name, 'filterEpoch');

%----------------   ICA START _-----------------------

fprintf('\nSTAGE 3: ICA filtering\n\n');
%filtered_epoch_eeg = pop_importdata('data', filtered_epoch);
%display(size(filtered_epoch_eeg.data));

ICA_input=pop_reref(epochEEG,[], 'keepref', 'on'); % Reference the EEG data to the average


%eegplot(ICA_input.data,'eloc_file',electrode_filename);
%-----------------------------------------------------------

% PCA reduce components

nchans=size(ICA_input.data,1);
ntimes=size(ICA_input.data,2);
nepochs=size(ICA_input.data,3);

ICA_input.data = reshape(ICA_input.data, [nchans ntimes*nepochs]);

k=25;
num_of_components_pca = floor(sqrt((ICA_input.pnts/k)));
%num_of_components_pca = min(num_of_components_pca, ICA_input.nbchan);

[ica_data, A, W] = fastica(ICA_input.data);   

ICA_input.icaweights=W;
ICA_input.icaact=ica_data;
ICA_input.icawinv=pinv(W);
ICA_input.icachansind=1:CHANNELS;

%display(ICA_input);
icaEEG.icaweights=W;
icaEEG.icaact=ica_data;
icaEEG.icawinv=pinv(W);
icaEEG.icachansind=1:CHANNELS;

lopass_freq=50;
[EEGdata] = filter_components(ICA_input, ICA_input, ica_data,lopass_freq);

ICA_input.data = EEGdata;
x=1:size(EEGdata,2);
x=x-FRAMES;
figure(fig_count),plot(x,EEGdata'), axis([-FRAMES size(EEGdata,2)-FRAMES -15 15]);
fig_count=fig_count+1;
title('Stage 3:Filtered EEG data after channel ICA Stage');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

[FileName, PathName] = uiputfile('./saved_data_outputs/*.mat', 'Save filtered EEG data after ICA');
name = strcat(PathName,FileName);
save(name, 'EEGdata');

%chans = [2 ];
%plot(ICA_input.data(chans,:)');


%--------------------------------------------------------------
% CHANNEL EPOCH FILTERING METHOD : LAST STAGE
%--------------------------------------------------------------

fprintf('\nSTAGE 4: SINGLE CHANNEL EPOCH FILTERING\n\n');
FRAMES=200;

epoched_data = epoch_data(FRAMES, TRIALS, EEGdata);
epochEEG.data = epoched_data;
epochEEG.pnts = size(epoched_data,2);
epochEEG.srate = SRATE;
epochEEG.trials = size(epoched_data,3);

for v=1:size(epochEEG.data,3)
    lengths_ep{v} = channel_epoch_filter(epochEEG,v);
    %display(v);
    %display(lengths_ep{v});
    %lengths_ep{v}=eeg_chans(logical(min_z(list_properties)));
end

ext_chans=[];
channels_epoch_eeg=h_epoch_interp_spl(epochEEG,lengths_ep,ext_chans);    % THIS IS THE FINAL EEG OUTPUT ON FASTER AFTER SINGLE CHANNEL EPOCHS
%channels_epoch_eeg = channel_epoch_filter(epochEEG);
final_data =epoch2eeg(channels_epoch_eeg.data);

%--------------------------------------------------
% Baseline correction
%--------------------------------------------------

fprintf('\n Performing baseline corrections\n');
data_average_baseline = mean(final_data(:,1:FRAMES));
data_correct = repmat(data_average_baseline, CHANNELS, 1);

%disp(size(data_correct));
for i=1:TRIALS
    final_data(:, (i-1)*FRAMES+1:i*FRAMES) = final_data(:, (i-1)*FRAMES+1:i*FRAMES) - data_correct;
end

x=1:size(final_data,2);
x=x-FRAMES;
figure(fig_count),plot(x,final_data'), axis([-FRAMES size(final_data,2)-FRAMES -15 15]);
fig_count=fig_count+1;
title('Stage 4:Filtered EEG data after baseline correction');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label

[FileName, PathName] = uiputfile('./saved_data_outputs/*.mat', 'Save final filtered data');
name = strcat(PathName,FileName);
save(name, 'final_data');

%{
%EEG = channels_epoch_eeg; % THIS IS THE FINAL EEG OUTPUT AFTER SINGLE CHANNEL EPOCHS
%display(channels_epoch_eeg);

%}
%}
%}

end
