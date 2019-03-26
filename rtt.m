function RTT = rtt(ms_array)
%% Calculate the Ratio of total time covered (RTT)
%  
%  Inputs : 1) ms_array   :EEG Microstate array of 1 epoch(a*b vector)
%                          a - number of epoches   b -  Number of time frame
%  Outputs: 1) RTT        :Percent of time covered by one template map within a given
%                          EEG epoch(1*a vector)
%                          a -  Number of microstates
%
%  Code: Soupee Li, Allard Shi, Xian Jiaotong University

%% Initialize the info. and Validation
[N_epoch, N_frame] = size(ms_array);
N_ms = length(unique(ms_array))-1;          
label = 1:N_ms;
RTT = zeros(N_epoch,N_ms);

%% Calculate RTT
for i = 1:N_ms
    for j = 1:N_epoch
        RTT(j,i) = length(find(ms_array(j,:) == label(i)))/(N_frame - length(find(ms_array(j,:) == 0)));
    end
end

end
