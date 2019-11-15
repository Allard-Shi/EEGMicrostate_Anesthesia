function ms_visualize(Ms,channel_loc)
%% Visualize micorstates
% 
%  Inputs: 1) Ms: Microstate matrix. (M*N matrix)
%              M - number of electrodes 
%              N - number of microstates
%          2) channel_loc: Location of channels (EEGLAB format)
%
% Author: Allard Shi, Xian Jiaotong University

%% Plot EEG microstates
[~,num_ms] = size(Ms);
for i =1:num_ms
    subplot(1,num_ms,i)
    topoplot(Ms(:,i)',channel_loc);
end
