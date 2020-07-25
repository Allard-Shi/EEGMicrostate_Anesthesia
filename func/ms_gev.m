function GEV = ms_gev(eegdata, ms_array, ms, choose_gfp)
%% Compute global explained varience across EEG signals
%
%  Inputs:  1) eegdata    :EEG data array (a*b*c matrix)
%                          a - Number of electrodes  b - Number of time frame 
%                          c - Number of epoches
%           2) ms_array   :Labeled Microstate EEG array (a*b matrix)
%                          a -  Number of time frame  b - Number of epoches
%           3) ms         :Microstates (a*b matrix)
%                          a - Number of electrodes  b - Number of microstates
%           4) choose_gfp :Choose GFP or not
%                         0: whole EEG series    otherwise: GFP peaks 
%
%  Outputs: 1) GEV        :Global field power within each microstate (a*1 vector)
%                         a -  Number of microstates
%
%  Anthor: Allard Shi  

%% Initialize the params
[~,N_frame, N_epoch] = size(eegdata);
N_ms = size(ms,2);
S = 0;
S_ms = zeros(N_ms,1);

%% Compute GEV
if choose_gfp == 0
    for i = 1:N_epoch
        for j = 1:N_frame
            if ms_array(j,i) == 0 
                continue
            else
               [C,c] = max(compute_spatial_correlation(eegdata(:,j,i),ms));
               S = S + sum(eegdata(:,j,i).*eegdata(:,j,i));       
               S_ms(c) = S_ms(c) + C^2 * sum(eegdata(:,j,i).*eegdata(:,j,i));
            end
        end
    end
else
    GFP = std(eegdata,1); 
    for i = 1:N_epoch
        [~,locs] = findpeaks(double(GFP(:,:,i)));
        for j = 1:length(locs)
            if ms_array(locs(j),i) == 0
                continue
            else
                C = compute_spatial_correlation(eegdata(:,locs(j),i),ms(:,ms_array(locs(j),i)));
                S = S + sum(eegdata(:,locs(j),i).*eegdata(:,locs(j),i));
                S_ms(ms_array(locs(j),i)) = S_ms(ms_array(locs(j),i)) + C*C * sum(eegdata(:,locs(j),i).*eegdata(:,locs(j),i));
            end
        end
    end
end
GEV = S_ms/S;

end