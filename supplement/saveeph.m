function saveeph(savefilename,thedata,varargin)
% saveeph: saves data as a Cartool evoked potential data file (.ep(h/sd/se))
%
% inputs: full path and name of the file to save, data as a 2-D numeric
% array where dimension 1 contains the timeframes, dimension 2 contains the
% channels, (optional) samplingrate as 1-D numeric array passed as
% 3rd argument
%
% outputs: .ep(h/sd/se) file
%
% Cartool: http://brainmapping.unige.ch/Cartool.htm
%
% author of this script: pierre.megevand@medecine.unige.ch


if strcmp(savefilename(end-3:end),'.eph')==1
    numtimeframes=size(thedata,1);
    numchannels=size(thedata,2);
    samplingrate=varargin{1};
    theheader=[numchannels numtimeframes samplingrate];
    dlmwrite(savefilename,theheader,'\t');
    dlmwrite(savefilename,thedata,'delimiter','\t','-append');
elseif strcmp(savefilename(end-2:end),'.ep')==1
    dlmwrite(savefilename,thedata,'delimiter','\t','-append');
elseif strcmp(savefilename(end-4:end),'.epsd')==1|strcmp(savefilename(end-4:end),'.epse')==1
    if nargin==2
        dlmwrite(savefilename,thedata,'delimiter','\t','-append');
    elseif nargin==3
        numtimeframes=size(thedata,1);
        numchannels=size(thedata,2);
        samplingrate=varargin{1};
        theheader=[numchannels numtimeframes samplingrate];
        dlmwrite(savefilename,theheader,'\t');
        dlmwrite(savefilename,thedata,'delimiter','\t','-append');
    end
end