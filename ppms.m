function Ppms = ppms(eeg, ms_array, N_ms, Fs)
%% Calculate Microstate-based Power
%  
%  Inputs : 1) eeg        : EEG time frame (a*b*c matrix)
%                           a - Number of electrodes
%                           b - Number of time frame
%                           c - Number of epoches
%           2) ms_array   : EEG microstate array of one epoch(a*b matrix)
%                           a - number of epoches  
%                           b - Number of time frame
%           3) N_ms       : Number of microstates
%           4) Fs         : Sampling rate (Hz)
%
%  Outputs: 1) Ppms       : microstate power(1*a vector)        
%                           a - Number of microstates
%
%  Code: Allard Shi, Xian Jiaotong University

%% Init params
Ppms  = zeros(1,N_ms);
[~,N_time,N_epoch] = size(eeg);
L_ms = zeros(1,N_ms);

%% Calculate Ppms
for i = 1:N_epoch
    for j = 1:N_time
        if ms_array(j,i) == 0
            continue
        else
            L_ms(1,ms_array(j,i)) = L_ms(1,ms_array(j,i)) + 1;
            Ppms(1,ms_array(j,i)) = Ppms(1,ms_array(j,i)) + sum(eeg(:,j,i).^2);
        end 
    end
end

Ppms = Ppms/sum(L_ms);

end