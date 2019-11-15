function MMD = mmd(ms_array, ms, Fs)
%% Compute the Mean Microstate Duration (MMD)
%  
% Inputs:  1)ms_array    :Microstate labeled array
%          2)ms          :Microstates
%          3)Fs          :Sampling frequency 
%
% Outputs: 1)MMD         :Mean microstate duration (ms)
%
% Anthor: Soupee Li, Allard Wen Shi  
% Copyright (C) Soupee Li, Allard Shi, Xian Jiaotong University 

%% Initialize the params
[N_num ,N_epoch] = size(ms_array);
N_ms = size(ms,2);                           % Number of microstates
ms_label = 1:N_ms;
MMD = zeros(N_epoch,N_ms);               
N_label = zeros(N_epoch,N_ms);               % Sum points of each microstate in every epoch     
N_period = zeros(N_epoch,N_ms);              % Nunber of microstate sections

%% Compute MMD
for i = 1:N_epoch                            
    for j = 1:N_ms                        
        N_label(i,j) = length(find(ms_array(:,i) == ms_label(j)));
        for k = 1:N_num
            if k > 1                        
                if ms_array(k,i) == ms_label(j) && ms_array(k-1,i) ~= ms_label(j)
                    N_period(i,j) = N_period(i,j)+1;
                end
            else
                if ms_array(k,i) == ms_label(j)
                    N_period(i,j) = N_period(i,j)+1;
                end
            end
        end
    end
end
MMD = N_label./N_period/Fs*1000;
 
end
