function [band_array, Corr, band_Cd] = bandicap(eeg_freq, ms_array, ms, threshold)
%% Fit microstates energy into different frequency bands
% 
% Inputs  £º1) eeg_freq  : EEG band signals (a*b cell)
%                          a - number of band 
%                          b - two col. 
%           2) ms_array  : microstates array (a*b matrix)
%                          a - number of time frame
%                          b - number of epoch
%           3) ms        : Number of microstates (a*b matrix)
%                          a - number of electrodes 
%                          b - number of microstates
%           4) threshold : Correlation threshold (dismiss Corr < threshold)
%                          
% Outputs:  1) band_array: fitting microstate to different bands (a*b matrix)
%                          a - number of time series 
%                          b - number of epoches
%           2) Corr      : Correlation (a*b matrix)
%                          a - number of time series 
%                          b - number of epoches
%
%           3) band_Cd   : microstate energy in different frequency bands (a*b matrix)
%                          a - number of band
%                          b - number of microstates
%
%
%  Anthor: Allard Wen Shi  
%  Copyright (C) EEGLab, Allard Shi, Xian Jiaotong University

%% Initialize the params
[N_band,~] = size(eeg_freq);
[N_p, N_epoch] = size(ms_array);
[~, N_ms] = size(ms); 
band_array = zeros(N_p, N_epoch);
Corr = zeros(N_p, N_epoch);
band_Cd = zeros(N_band, N_ms); 

%% Compute band info. 
for i = 1:N_epoch
    C =zeros(N_band, N_p);
    for j = 1:N_band
        for k = 1:N_p
            if ms_array(k,i) == 0
                C(j,k) = 0;
                continue;
            end
            C(j,k) = compute_spatial_correlation(eeg_freq{j,2}(:,k,i),ms(:,ms_array(k,i)));   
        end
    end
    [Corr(:,i),band_array(:,i)] = max(C,[],1);
    Corr(find(Corr(:,i) < threshold),i) = -1;
    band_array(find(Corr(:,i) < threshold),i) = -1;
end

for i = 1:N_band
    for j = 1: N_ms
        if length(find(ms_array == j & Corr~= 0 & Corr~=-1)) == 0
            band_Cd(i,j) = 0;
        else
            band_Cd(i,j) = length(find(band_array == i & ms_array ==j))/length(find(ms_array == j & Corr~= 0));
        end
    end
end

band_Cd = reshape(band_Cd',1,N_band*N_ms);
end