function [ms_array, Corr] = ms_fit(eegdata, ms, lamda, b, varargin)
%% Fitting eegdata raw data using microstates
%  
%  Inputs:  1) eegdata   : EEG data array (a*b*c matrix)
%                          a - Number of electrodes  b - Number of time frame 
%                          c - Number of epoches
%           2) ms        : Microstates (a*b matrix)
%                          a - Number of electrodes  b - Number of microstates
%           3) lamda     : Penalty parameter for nonsmoothness
%           4) b         : Window semi-length
%           5) default   : tol (1e-6) 
%                          corr_thrsh - correaltion threshold (0.00 - 1.00) 
%
%  Outputs: 1) ms_array  :Labeled Microstate EEG array (a*b matrix)
%                         a -  Number of time frame  b - Number of epoches
%           2) Corr      :Maximum of spatial similarity between microstate
%                         topography and EEG topography (1*b cell)
%                         b - Number of epoches
%
%  Anthor: Soupee Li, Allard Wen Shi  

%% Default Params
if ~exist('tol'),
    tol = 1e-6;
end
if ~exist('corr_thrsh'),
    corr_thrsh = 0.50;
end

%% Initialize the info. and Validation 
GFP = std(eegdata,1);                               % Recalculate GFP values
[N_e, N_ms] = size(ms);                             % Number of electrodes and microstates
[~,N_num,N_epoch] = size(eegdata);                  % Time frame and Number of epoches                      

if size(eegdata) ~= N_e
    fprintf('Number of electrodes in micorstates and EEG does not match!\n');
end

%% Find peak locations in eegdata and compute correlation of each microstate.
for i = 1:N_epoch            
    [~,locs] = findpeaks(double(GFP(:,:,i)));       % Find the peaks      
    Locs{i} = locs';    
    [Corr{i},C{i}] = max(compute_spatial_correlation(eegdata(:,locs,i),ms),[],2);
end
Ka = C;

%% Smoothing Algorithms
for n = 1:N_epoch
    sigma_o = 0; sigma_u = 1000;
    frame = length(Locs{n});
    Nbkt = zeros(frame-2*b,N_ms);                   % Define window
    Peak_map = eegdata(:,Locs{n},n);
    e = 0;
    for j = 1:frame
        e = e + (Peak_map(:,j)'*Peak_map(:,j)-(ms(:,C{n}(j))'*Peak_map(:,j))^2);     
    end
    e = e/(frame*(N_e-1));
    
    while abs(sigma_o-sigma_u) > tol*sigma_u
        sigma_o = sigma_u;
        
        % Compute Nbkt between a window whose length is 2b+1
        for t = 1+b:frame-b
            for k = 1:N_ms
                Nbkt(t-b,k) = length(find(C{n}(t-b:t+b) == k));
            end

            % argmin(k)
            S = zeros(N_ms,1);
            for k = 1:N_ms
                S(k) = (Peak_map(:,t)'*Peak_map(:,t)-(ms(:,k)'*Peak_map(:,t))^2)/(2*e*(N_e-1))-lamda*Nbkt(t-b,k);
            end
            [~,Min_index] = min(S);
            Ka{n}(t) = Min_index; 
        end
        C{n} = Ka{n};

        % Compute new sigma_u
        sigma_u = 0;
        for j = 1:frame
            sigma_u = sigma_u + (Peak_map(:,j)'*Peak_map(:,j)-(ms(:,C{n}(j))'*Peak_map(:,j))^2);     
        end
        sigma_u = sigma_u/(frame*(N_e-1));
    end
    
end

%% Label process

% some studies will set thresholds for correlation. 
% for i = 1:N_epoch
%     Loc_lowcorr = find(Corr{i} < corr_thrsh);
%     C{i}(Loc_lowcorr) = 0;
% end

ms_array=zeros(N_num,N_epoch);                       % Initialize the labeled array
for i = 1:N_epoch
    inner_index = 1;
    for j = 1:N_num
        if j < floor(1 + Locs{i}(2)/2)               % Discard the first and last microstate arrays
            ms_array(j,i) = 0;
        elseif j > floor((N_num+Locs{i}(end -1))/2)
            ms_array(j,i) = 0;
        elseif inner_index < size(Locs{i},1) && j == floor((Locs{i}(inner_index) + Locs{i}(inner_index+1))/2)
            ms_array(j,i) = C{i}(inner_index+1);
            inner_index = inner_index + 1;
        else
            ms_array(j,i) = C{i}(inner_index);
        end
    end 
end

end