function [eeg_data] = Makinen1a(channels, epochs, time_points)

% function Makinen1a
%
% function plots the Figure 1a from the Makinen et al. (2005)

trials = epochs;
frames = time_points;
chans = channels;

for i=1:chans
    mysig = Makinen (frames, trials, 1000, 50);
    mysig = reshape (mysig, frames, trials);
    k(i,:) = mean(mysig');
end

eeg_data = k;
disp(size(k));
plot (k');
ylabel ('ERP');

%mysig = reshape (mysig, 400, trials);
%{
subplot (3,1,1);
plot (mysig');
disp(size(mysig))
ylabel ('EEG');
subplot (3,1,2);
plot (mean(mysig'));
ylabel ('ERP');
subplot (3,1,3);
plot (var(mysig'));
ylabel ('Variance');
%}