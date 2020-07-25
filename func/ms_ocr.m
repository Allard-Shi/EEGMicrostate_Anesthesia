function ocr_v = ms_ocr(ms_array, N_ms, Fs, window)
%% Calculate the Occurrence of specific microstates
%  
%  Inputs:  1) ms_array: Microstate array(a*b matrix)
%                         a - Number of time frame
%                         b - Number of epochs
%           2) N_ms    : Number of microstates
%           3) Fs:       Sampling rate (Hz)
%           4) window  : Length of window (sec)
%
%  Outputs: 1) ocr_v:    Occurrence, the mean number of distinct microstates of a
%                        given class occuring within a 1 sec window (a*1 vector)
%                         a - Number of microstates
%
%  Anthor: Soupee Li, Allard Wen Shi  

%% Initialize the info.  
[N_frame, N_epoch] = size(ms_array);
ocr_v = zeros(N_ms,1);
N_num = sum(sum(ms_array ~= 0));                  % Aquire the length of ms array without zero

%% Calculate the Occurrence
for i = 1:N_epoch
    temp =  ms_array(1,i);
    for j = 2:N_frame
        if temp ~= ms_array(j,i) && ms_array(j,i)~=0
            ocr_v(ms_array(j,i),1) = ocr_v(ms_array(j,i),1) + 1;
        end
        temp = ms_array(j,i);  
    end
end
ocr_v = ocr_v*window*Fs/N_num;

end