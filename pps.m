function PPS = pps(eeg,Fs)
%% Global field power peaks per second (PPS)
% 
% Inputs  £º1) eeg     : EEG signals (a*b*c matrix)
%                        a - number of electrodes  
%                        b - number of time series  
%                        c - number of epoches
%           2) Fs      : Sampling frequency (Hz)
%
% Outputs:  1) PPS     :Global field power peaks per second (times/second)
%
% Code: Allard Shi, Xian Jiaotong University

%% Initialize the params
GFP = std(eeg,1);
[~,N_p,N_epoch] = size(GFP);
L = N_p/Fs;
pps_all = zeros(1,N_epoch);

%% Compute PPS
for i = 1:N_epoch
     [~,locs] = findpeaks(double(GFP(:,:,i)));
     pps_all(1,i) = length(locs)/L;
end
PPS = mean(pps_all);

end