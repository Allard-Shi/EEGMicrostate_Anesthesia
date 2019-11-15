%% Fit Cartool
clear;clc;
filepath = '';
resultpath = '';
datainfo_name = '\datainfo.mat';
sef_savepath = '';
load(strcat(filepath,datainfo_name));
N_people = 20;                                              % Number of People

%% create electrode file for CarTool
for n = 1:N_people
    % obtain the name of file in the same section
    filename = strcat(datainfo(4*n,1),'.set'); 
    eeg = pop_loadset(filename, filepath);
    fprintf('Import People %d  Success!\n',n);
    EEG = eeg.data;
    Fs = eeg.srate; T = 1/Fs;                 % sampling frequency
    N_p = eeg.pnts;                           % Number of sampling points
    N_e = eeg.nbchan;                         % Number of electrodes
    N_trial = eeg.trials;                     % Number of trials
    E_Loc = eeg.chanlocs;                     % Location of electrodes
    samplingrate = Fs*ones(N_e,1);            % resampling rate matrix
    for i = 1:N_trial
        savesef([sef_savepath,num2str(datainfo{4*n-3,2}),'\',num2str(n),'\',num2str(i),'.sef'],EEG(:,:,i)',samplingrate); 
    end  
end 

for i = 1:length(E_Loc)
    loc{i,1} = num2str(E_Loc(i).X,6);
    loc{i,2} = num2str(E_Loc(i).Y,6);
    loc{i,3} = num2str(E_Loc(i).Z,6);
    loc{i,4} = E_Loc(i).labels;
end

% Save data in EEG_locs.els;
fidpath = 'EEG_locs.els';
fid = fopen(fidpath,'wt');
fprintf(fid,['ES01\n' int2str(size(loc,1)) '\n1\n10-10 System\n' int2str(size(loc,1)) '\n3\n']);
for i = 1 : size(loc,1)
    fprintf(fid,'%s\t %s\t %s\t %s\t\n',loc{i,:});
end
fclose(fid);
fprintf('Electrode file saved!\n');