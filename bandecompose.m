function [Band_power,Band_ms,Spectrum] = bandecompose(memdimf,ms_array,band,N_ms,Fs)
%% Build up relationship between microstates and frequency band 
%
% Inputs  £º1) memdimf    - Intrinsic mode function (IMFs) after MEMD (1*a cell)
%                           a - number of epoch
%           2) ms_array   - microstates array (a*b matrix)
%                           a - number of time frame
%                           b - number of epoch
%           3) band       - band information (a*b matrix)
%                           a - number of selected band
%                           b - first number  - low  freq. 
%                               second number - high freq.
%           4) N_ms       - Number of microstates 
%           5) Fs         - Sampling rate
%                          
% Outputs:  1) Band_power - Signal power distributed to different band (a*b*c matrix)
%                           a - number of time series 
%                           b - number of epoches
%                           c - number of selected band
%           2) Band_ms    - Band power distributed to different microstates (a*b matrix)
%                           a - number of selected band
%                           b - number of microstates
%           3) Spectrum   - Power spectrum: 
%
%  Anthor: Allard Wen Shi  
%  Copyright (C) Allard Shi, Xian Jiaotong University

%% init params
N_epoch = size(memdimf,2);
[N_band,~] = size(band);
[N_e,~,N_time] = size(memdimf{1}); 
Band_energy = zeros(N_time,N_epoch,N_band);
Band_ms = zeros(N_band,N_ms);

%% HHT analysis
for i = 1:N_epoch
    imf = memdimf{1,i};
    for j = 1:N_e
        if j == 1
            [P,~,~,imfinsf,imfinse] = hht(permute(imf(j,:,:),[3,2,1]),Fs);
            P = full(P);
        else
            [P1,~,~,imfinsf,imfinse] = hht(permute(imf(j,:,:),[3,2,1]),Fs);
            P = P + full(P1);
        end
        %[~,~,~,imfinsf,imfinse] = hht(imf,Fs);
        N_imf = size(imfinsf,2);
        for k = 1:N_time
            for l = 1:N_imf
                for m = 1:N_band
                    if imfinsf(k,l) >= band(m,1) && imfinsf(k,l) < band(m,2)
                        Band_energy(k,i,m) = Band_energy(k,i,m) + imfinse(k,l);
                    end
                end
            end
        end
    end
    Spectrum{1,i} = P;
end

%% Transfer to microstates
temp = 0;
L_ms = zeros(N_ms,1);
for i = 1:N_epoch
    for j = 1:N_time
        if ms_array(j,i) == 0
            temp = temp + 1;
            continue
        else
            for k = 1:N_band
                Band_ms(k,ms_array(j,i)) = Band_ms(k,ms_array(j,i)) + Band_energy(j,i,k);
            end        
        end
    end
    for j = 1:N_ms
        L_ms(j) = L_ms(j) + length(find(ms_array(:,i) == j));
    end
end

for i = 1:N_ms
    for j = 1:N_band
        Band_ms(j,i) = Band_ms(j,i)/(sum(L_ms));
    end
end

end