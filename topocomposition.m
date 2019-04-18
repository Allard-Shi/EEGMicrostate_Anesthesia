function [ms_energy,ms_power] = topocomposition(memdimf,ms_array,band,N_ms,Fs)
%% Build up relationship between microstate topography and frequency band 
%
% Inputs  £º1) memdimf   : Intrinsic mode function (IMFs) after MEMD (1*a cell)
%                          a - number of epoch
%           2) ms_array  : microstates array (a*b matrix)
%                          a - number of time frame
%                          b - number of epoch
%           3) band      : band information (a*b matrix)
%                          a - number of selected band
%                          b - [low  freq.  high freq.]                
%           4) N_ms      : Number of microstates 
%           5) Fs        : Sampling rate
%                          
% Outputs:  1) ms_energy : Signal energy distributed to different band and topographies (a*b*c matrix)
%                          a - number of electrodes 
%                          b - number of bands
%                          c - number of microstates
%           2) ms_power  : Band power distributed to different band and topographies (a*b*c matrix)
%                          a - number of selected band
%                          b - number of microstates
%                          c - number of microstates
%
%  Anthor: Allard Wen Shi  
%  Copyright (C) Allard Shi, Xian Jiaotong University

%% Initialize the params
N_epoch = size(memdimf,2);
[N_band,~] = size(band);
[N_e,~,N_time] = size(memdimf{1}); 
ms_energy = zeros(N_e,N_band,N_ms);
L = 0;

%% MEMD-HT analysis
for i = 1:N_epoch
    imf = memdimf{1,i};
    for j = 1:N_e
        [~,~,~,imfinsf,imfinse] = hht(permute(imf(j,:,:),[3,2,1]),Fs);
        N_imf = size(imfinsf,2);
        for k = 1:N_time
            if ms_array(k,i) == 0
                continue
            else
                L = L+1;
                for m = 1:N_band
                    for n = 1:N_imf
                        if imfinsf(k,n) >= band(m,1) && imfinsf(k,n) < band(m,2)
                             ms_energy(j,m,ms_array(k,i)) = ms_energy(j,m,ms_array(k,i)) + imfinse(k,n);
                        end 
                    end
                end
            end
        end
    end
end
ms_power = ms_energy / L;

end