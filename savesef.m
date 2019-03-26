function savesef(savefilename,thedata,samplingrate)
% savesef: saves data as a Cartool simple EEG data file (.sef)
%
% inputs: full path and name of the file to save, data as a 2-D numeric
% array where dimension 1 contains the timeframes, dimension 2 contains the
% channels, samplingrate as 1-D numeric array
%
% outputs: .sef file
%
% Cartool: http://brainmapping.unige.ch/Cartool.htm
%
% author of this script: pierre.megevand@medecine.unige.ch


% define fixed part of header
version='SE01';
numchannels=size(thedata,2);
numauxchannels=0;
numtimeframes=size(thedata,1);
year=0;
month=0;
day=0;
hour=0;
minute=0;
second=0;
millisecond=0;

% open savefilename for writing
fid=fopen(savefilename,'w');

%write fixed part of header
fwrite(fid,version,'int8');
fwrite(fid,numchannels,'int32');
fwrite(fid,numauxchannels,'int32');
fwrite(fid,numtimeframes,'int32');
fwrite(fid,samplingrate,'float32');
fwrite(fid,year,'int16');
fwrite(fid,month,'int16');
fwrite(fid,day,'int16');
fwrite(fid,hour,'int16');
fwrite(fid,minute,'int16');
fwrite(fid,second,'int16');
fwrite(fid,millisecond,'int16');

% define and write variable part of header
for i=1:numchannels
    fwrite(fid,101,'int8');
    currentchannel=uint8(num2str(i));
    for j=1:size(currentchannel,2)
        fwrite(fid,currentchannel(j),'int8');
    end
    for k=j+2:8
        fwrite(fid,0,'int8');
    end
end

% write data
count = fwrite(fid,thedata','float32');

% close file
fclose(fid);