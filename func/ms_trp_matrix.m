function TP = ms_trp_matrix(gfp,ms_array,ms)
%% Calculate transition probability matrix
%  
%  Inputs:  1) gfp       : Global field power (1*b*c matrix)
%                          b - Number of time frame c - Number of epoches
%           2) ms_array  : Microstate frame (a*b matrix)
%                          a - Number of time frame  b - Number of epoches
%           3) ms        : Microstates (a*b matrix)
%                          a - Number of electrodes  b - Number of microstates
%
%  Outputs: 1) TP        : The transition probability matrix of input epoch([a*a,c] matrix)
%                          a - Number of microstate 
%                          c - Number of Epoches
%
% Anthor: Soupee Li, Allard Shi  

%% Initialize the info. and Validation
[~, N_ms] = size(ms);                                          % Number of electrodes and microstates
[~, N_epoch] = size(ms_array);                                 % Time frame and Number of epoches                      
tp = zeros(N_ms, N_ms);


%% Calcualte the transition Probability Matrix
for i = 1:N_epoch
    [~,locs] = findpeaks(double(gfp(:,:,i)));                  % Find the peaks      
    array = ms_array(locs,i);
    array(array == 0) = []; 
    temp = array(1);
    for j = 2:length(array)-1
        tp(temp,array(j)) =tp(temp,array(j)) + 1;
        temp = array(j);
    end  
end
tp_s = sum(tp,2);

% Compute final matrix
for j = 1:N_ms
    if tp_s(j) == 0
        continue
    else
        tp(j,:) = tp(j,:)/tp_s(j);
    end
end
TP = reshape(tp',N_ms*N_ms,1);
end