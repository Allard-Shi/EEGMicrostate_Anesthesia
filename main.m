%% Microstates Calculation and Validation
%
%  This file is aimed to compute the GFP and microstates of EEG 
%  The data used in the research has been avaliable in the following website.
%  link: https://www.repository.cam.ac.uk/handle/1810/252736
%
%  Anthor: Allard Shi, Xian Jiaotong University 

%% Main code
clear;clc;
addpath('func');
filepath = '';
resultpath = ''; 
datainfo_name = '';
load(strcat(filepath,datainfo_name));
N_people = 20;                                             % Number of People
Num_cluster = ; Loop = 100;                                % Set the maximum of the number of clusters
R_square = zeros(N_people,Num_cluster);
R_cv = zeros(Num_cluster,1);
W = zeros(Num_cluster,1);
MS = [];
b = 1; lamda = 0.05;                                           % Set smoothing params (smooth the microstate sequence)

%% Load data and compute the EEG microstates
iter = 1;
for m = 1:iter
    for c =
        for n = 1:N_people   
		
            % obtain the name of file in the same section
            filename = strcat(datainfo(4*n-3,1),'.set'); 
            eeg = pop_loadset(filename, filepath);
            fprintf('Import People %d  Success!\n', n);
            EEG = eeg.data;
            Fs = eeg.srate; T = 1/Fs;                           % sampling frequency
            N_p = eeg.pnts;                                     % Number of sampling points
            N_e = eeg.nbchan;                                   % Number of electrodes
            N_trial = eeg.trials;                               % Number of trials
            E_Loc = eeg.chanlocs;                               % Location of electrodes  
            GFP = std(EEG,1);                                   % Define Global field power (preprocessed mean = 0)
            Kappa = [];

            % Extract all GFP peak maps
            Gamma{n} = zeros(N_trial, N_e, 4);
            MS = [];
            for i = 1:N_trial
                [X,locs] = findpeaks(double(GFP(:,:,i)));
                [Kappa_max,gamma,R_temp,R_cv_temp,W_temp,M_temp] = mkmeans95(EEG(:,locs,i), ,Loop,b,lamda);
                Kappa{i} = Kappa_max;
                MS = [MS,gamma];
            end
            K{n} = Kappa;
            fprintf('People %d Finished!\n', n);
        end
		
		% Second clustering
		[L{c},G{c},R{c},R_cv(c,m),W(c,m),~] = mkmeans95(MS, ,Loop,b,lamda);
    end
end
