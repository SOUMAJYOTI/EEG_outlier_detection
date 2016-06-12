function [epoched_data] = epoch_data(FRAMES, TRIALS, data)

% Inputs:
%  FRAMES - number of signal frames per each trial
%  TRIALS - number of simulated trials
%  SRATE - sampling rate of simulated signal
% Output:
%  synthetic_data - simulated EEG data; vector: 32 by frames*epochs containing concatenated trials

%%

SRATE = 250;

for i=1:TRIALS
	epoched_data(:,:,i) = data(:,(i-1)*FRAMES+1:i*FRAMES);
end

end

