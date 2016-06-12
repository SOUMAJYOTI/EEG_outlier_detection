function test_faster(TRIALS, CHANNELS, NAME_DATA)

% test_faster() - Processes all the stages in the EEG artifact rejection
%
% Usage:
%  >> test_faster(5, 64, 'noisy_datas/1.mat'); 
%
% Inputs:
%   TRIALS- number of epochs
%   CHANNELS - number of channels

%clear all;
clc
SRATE =250;
fig_count=1;
FRAMES = 200;

TRIALS=TRIALS+1;
load(NAME_DATA);

original_signal(:,:) = save_data(1,:,:);
noisy_signal(:,:) = save_data(2,:,:);

disp(size(original_signal));
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

EEG = pop_importdata('data',noisy_signal);
%disp(EEG)
if CHANNELS==32
    electrode_filename='data_sample/32.ced';
elseif CHANNELS==64
    electrode_filename='data_sample/64.ced';
else
    electrode_filename='data_sample/128.ced';
end

EEG.chanlocs=readlocs(electrode_filename);
EEGdata=EEG.data;

%rawEEG= EEG;

eeg_chans = 1:EEG.nbchan;
ref_chan = [];
ext_chans=[];

EEG = pop_reref( EEG, [], 'keepref', 'on');  % Average referenced
epochEEG=EEG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpf_on=0;
lpf_on=0;
notch_on=0;

do_hipass=hpf_on;
do_lopass=lpf_on;
do_notch=notch_on;

hpf_freq=1;
hpf_bandwidth=0.50;
hpf_ripple=0.05;
hpf_attenuation=80;

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
%------------------------
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

% Notch filter
%---------------------
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
%THEIR FASTER METHOD FOR CHANNEL FILTERING
%-------------------------------------------------------------
% Channel interpolation options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 Automatic interpolation of bad channels on or off (1 / 0)
% 2 Radius for channel interpolation hypersphere (integer > 0)
% 3 Automatic interpolation of channels per single epoch at end of process (1 / 0)
% 4 Radius for epoch interpolation hypersphere (integer > 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('FASTER METHOD FOR CHANNEL FILTERING:');
    
eeg_chans = 1:tempEEG.nbchan;
ref_chan = [];

list_properties = channel_properties(tempEEG,eeg_chans,ref_chan);
lengths = min_z(list_properties);
chans_to_interp = eeg_chans(logical(lengths));
chans_to_interp = setdiff(chans_to_interp,ref_chan); % Ref chan may appear bad, but we shouldn't interpolate it!
display(chans_to_interp);

tempEEG = eeg_interp(tempEEG,chans_to_interp);
%name = strcat('saved_data_model/chanFiltered_faster_',num2str(noise_channel_perc),'_percent_noise_',num2str(noise_amp),'_noise_amp');
%save(name, 'tempEEG');
tempdata=tempEEG.data;
x=1:size(tempEEG.data,2);
x=x-FRAMES;
figure(fig_count),plot(x,EEG.data'), axis([-FRAMES (TRIALS-1)*FRAMES -15 15]);
title('Stage1: Channel filtered EEG');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;
[FileName, PathName] = uiputfile('./saved_data_outputs/channel_filtered_output/*.mat', 'Save channel filtered EEG data using FASTER');
name = strcat(PathName,FileName);
save(name, 'tempdata');
%fprintf('\nThe data output has been saved in the saved_data_model folder\n\n');
%disp(tempEEG);
EEG=tempEEG;

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
% EPOCH FILTERING METHOD USING FASTER 
%--------------------------------------------------------------
epdata = EEG.data;
epoched_data = epoch_data(FRAMES, TRIALS, epdata);

%disp(size(epoched_data));
epochEEG.data = epoched_data;
epochEEG.pnts = size(epoched_data,2);
epochEEG.srate = SRATE;
epochEEG.trials = size(epoched_data,3);

list_properties = epoch_properties(epochEEG,eeg_chans);
num_of_epochs = size(epochEEG.data, 3);
[lengths] = min_z(list_properties);

k=1;

epoch_noisy = [];
for i=1:num_of_epochs
    if lengths(i)==1
        epoch_noisy(k)=i;
        k=k+1;
    end
end

fprintf('The noisy epochs detected by FASTER are\n\n');
if ~isempty(epoch_noisy)
    display(epoch_noisy);
else
    fprintf('\tNo noisy epochs detected!!\n');
end
k=1;
for i=1:size(epochEEG.data,3)
    if ismember(i,epoch_noisy)==0
        filtered_epoch_faster(:,:,k) = epochEEG.data(:,:,i); 
        k=k+1;
    end
end

%disp(size(filtered_epoch_faster));
TRIALS=size(filtered_epoch_faster,3);
filtered_epoch=epoch2eeg(filtered_epoch_faster);
tempEEG.data = filtered_epoch;
EEG=tempEEG;
EEG.pnts=size(filtered_epoch_faster,2);
EEG.trials=1;
EEG=pop_reref(tempEEG,[], 'keepref', 'on'); % Reference the EEG data to the average
EEG.pnts=size(filtered_epoch,2);
%disp(EEG);

%fprintf('trials\n');
%disp(size(tempEEG.data));
%disp(size(filtered_epoch));

x=1:size(EEG.data,2);
x=x-FRAMES;
figure(fig_count),plot(x,EEG.data'), axis([-FRAMES size(EEG.data,2)-FRAMES -15 15]);
fig_count=fig_count+1;
title('Epoch filtered EEG');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;

[FileName, PathName] = uiputfile('./saved_data_outputs/epoch_filtered_output/.mat', 'Save epoch filtered EEG data using FASTER');
name = strcat(PathName,FileName);

%fprintf('trials\n');
%disp(size(filtered_epoch));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ICA Filtering stage
    % %%%%%%%%%%%%%%%%%%%
    % ICA options %
    %%%%%%%%%%%%%%%%%%%%%
    % 1 ICA on or off (1 / 0)
    % 2 Auto component rejection on or off (1 / 0)
    % 3 Radius for component rejection hypersphere (integer > 0)
    % 4 EOG channels (vector of integers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_ica=1;
k_value=25;
component_rejection_on=1;
keep_ICA=1;

do_ica = run_ica;
do_component_rejection = component_rejection_on;
EOG_chans = [];
ica_chans = 1:size(EEG,1);
 
%%%%%%%%%%
% Do ICA %
%%%%%%%%%%
 if do_ica && (keep_ICA || isempty(EEG.icaweights))
    %num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
    %num_pca = min(num_pca,length(setdiff(ica_chans,chans_to_interp)));
    k=25;
    num_pca = floor(sqrt((EEG.pnts/k)));
    disp(num_pca);
   % nchans=size(EEG.data,1);
    %ntimes=size(EEG.data,2);
    %nepochs=size(EEG.data,3);

    %EEG.data = reshape(EEG.data, [nchans ntimes*nepochs]);
    
    [weights sphere] = runica(EEG.data, 'pca', num_pca);
    EEG.icaweights = weights;
    EEG.icasphere = sphere;
    EEG.icawinv = pinv(weights);
    EEG.icachansind = 1:size(EEG.data,1);
    
 end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Component rejection section %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Also includes topoplots   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
lopass_freq=50;
list_properties = component_properties(EEG,EOG_chans,[lopass_freq-5 lopass_freq+5]);
%display(list_properties);
%o.ica_options.rejection_options.measure(2)=0;
        
[lengths] = min_z(list_properties);
bad_comps=find(lengths);
fprintf('\nThe bad component detected is: ');
disp(bad_comps);

%disp(size(EEG.data));
%disp(EEG);
% Reject
if ~isempty(find(lengths,1))
    EEG = pop_subcomp(EEG, find(lengths), 0);
end
%disp(size(EEG.data));
x=1:size(EEG.data,2);

figure(fig_count),plot(x,EEG.data'), axis([0 size(EEG.data,2) -15 15]);
title('ICA filtered EEG');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;


%--------------------------------------------------------------
% CHANNEL EPOCH FILTERING METHOD : LAST STAGE
%--------------------------------------------------------------
fprintf('\nPerforming filtering in single channel epochs\n\n');
eeg_chans=1:size(EEG.data,1);
ext_chans=[];

FRAMES=200;

%disp(EEG);
epoched_data = epoch_data(FRAMES, TRIALS, EEG.data);
epochEEG.data = epoched_data;
epochEEG.pnts = size(epoched_data,2);
epochEEG.srate = SRATE;
epochEEG.trials = size(epoched_data,3);
%disp(size(epochEEG.data,3));

for v=1:size(epochEEG.data,3)
    list_properties = single_epoch_channel_properties(epochEEG,v,eeg_chans);
    lengths_ep{v}=eeg_chans(logical(min_z(list_properties)));
end
%disp(size(epochEEG.data));
%disp(lengths_ep{2});
EEG=h_epoch_interp_spl(epochEEG,lengths_ep,ext_chans);    % THIS IS THE FINAL EEG OUTPUT ON FASTER AFTER SINGLE CHANNEL EPOCHS
EEG.data=epoch2eeg(EEG.data);

%----------------------------------------------------
%   BASELINE CORRECTION BY DEDUCTING THE AVERAGE
%----------------------------------------------------

fprintf('\nPerforming baseline correction\n');
data_average_baseline = mean(EEG.data(:,1:FRAMES));
data_correct = repmat(data_average_baseline, CHANNELS, 1);
for i=1:TRIALS
    EEG.data(:, (i-1)*FRAMES+1:i*FRAMES) = EEG.data(:, (i-1)*FRAMES+1:i*FRAMES) - data_correct;
end

x=1:size(EEG.data,2);
figure(fig_count),plot(x,EEG.data'), axis([0 size(EEG.data,2) -15 15]);
title('Stage 4: Final EEG after baseline correction');
xlabel('Time(ms)') % x-axis label
ylabel('Potential(uV)') % y-axis label
fig_count=fig_count+1;

fprintf('\n\n Completed all stages\n');
%}
%}
end
