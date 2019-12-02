function waveletSortMCDmultiple(fileName,saveFile)
%Reads a provided MCD file and returns the rising (positive) and falling
%(negative) edges of light (A1) and injection (A2) outputs. Also outputs
%the stimulus tags (D1), cycle tags (D2), pulse tags (D3), and time
%extracted from the MCD file

% MCD File components:
% A1 - Analog input 1 - Output from LED light source
% A2 - Analog input 2 - Output from multiport pressure injection system
% D1 - Digital input 1 - Stimulus tag indicating when each stimulus starts
% D2 - Digital input 2 - Cycle tag indicating the start of each stimulus trial
% D3 - Digital input 3 - Pulse tag indicating when each stimulus pulse
% begins for cases where multiple stimuli occur in one trial

spikeData = []; % Intiailize the spike data array
allDat = []; % initialize the data structure

for ii = 1:length(fileName) % Load each MCD file and extract data
    [cSpikeData,initialTime(ii),cDat] = readSpikeData(char(fileName(ii))); % Read spike data and initial time points from MCD file
    finalTime(ii) = cDat.finalTime; % Save the final time of each MCD file to reconcile time between files
    spikeData = appendSpikeData(spikeData,cSpikeData); % save spike data to single data structure
    allDat = appendStruct(allDat,cDat); % Save stimulus, cycle, pulse times, and analaog inputs into data structure
end


minTime = min(initialTime); % Find the starting time for the first MCD file
maxTime = max(finalTime)-minTime; % Find the final time of the last MCD file
for ii = 1:length(spikeData)
    spikeData(ii).timeStamps = spikeData(ii).timeStamps-minTime; % Correct all spike timestamps to the starting time for the first MCD file
end


A1 = allDat.A1; % Extract the rising/falling edges from the A1 input
A1(A1>0) = A1(A1>0)-minTime; % Correct rising edges times 
A1(A1<0) = A1(A1<0)+minTime; % Correct falling edges times

A2 = allDat.A2; % Extract the rising/falling edges from the A2 input
A2(A2>0) = A2(A2>0)-minTime; % Correct rising edge times
A2(A2<0) = A2(A2<0)+minTime; % Correct falling edges times

% A1 = allDat.A1-minTime;
% A2 = allDat.A2-minTime;
D1 = allDat.D1-minTime; % Correct all D1 input rising edges
D2 = allDat.D2-minTime; % Correct all D2 input rising edges
D3 = allDat.D3-minTime; % Correct all D3 input rising edges


%% Initialize Data Folder
% "folderName" should be set to a fast drive with at least 4 GB of space
% This folder will be used for the Wave_clust spike sorting


[ret, name] = system('hostname'); % Determine which computer is being utilized
% if strcmp(name(1:end-1),'TroyLabLenovo')
%     folderName = 'F:\Temp Folder';
% elseif strcmp(name(1:end-1),'Corey-PC')
%     folderName = loader('Processed Data/Temp Data Folder',1);
% elseif strncmp(name,'SaggereLab20643',15)
%     folderName = 'E:\temp';
% else
%     folderName = uigetdir;
% end

folderName ="D:\Troy\Wave_Clus UIC\20151111\Recording1";

% This code clears existing sorting files but is not usually necessary:

% d = dir(folderName);
% 
% cFileNames = arrayfun(@(x) d(x).name,1:length(d),'Unif',false)';
% cFileNames = cFileNames(~strcmp(cFileNames,'.')&~strcmp(cFileNames,'..'));
% for ii = 1:length(cFileNames)
%     delete(sprintf('%s/%s',char(folderName),char(cFileNames(ii))));
% end



oDir = cd; % Store original directory for future use

fileID = fopen(sprintf('%s/files.txt',folderName),'wt'); % Create a text file called "files.txt" that contains the names of all the spike data files for Wave_clus sorting

for ii = 1:length(spikeData) % Go through the spike data and save the spikes from each channel as a separate .mat file and include the filename in "files.txt" 
    cName = sprintf('%s/ch_%s.mat',char(folderName),char(spikeData(ii).channelName)); % Create the filename of the .mat file for the current channel of spikes
    index = spikeData(ii).timeStamps'; % extract the timestamps for the current channel
    spikes = spikeData(ii).waveforms'; % extract the spike waveforms for the current channel
    save(cName,'index','spikes'); % save the timestamps and spike waveforms to a .mat file
    if ~isempty(spikes) % if there are spikes in the current channel, write the filename into the "files.txt" text file
        cNameShort = sprintf('ch_%s.mat\n',char(spikeData(ii).channelName));
        fprintf(fileID,'%s',cNameShort); % write filename to "files.txt"
    end
end
fclose(fileID); % Close the "files.txt" file

%% Perform Sorting
close all
cd(folderName) % Navigate to the temporary file directory for spike sorting
Do_clustering(sprintf('%s/files.txt',folderName)); % Perform Wave_clust

close all
clc

%% Save info for wavelet clustering

for ii = 1:length(spikeData) % Go through the Wave_clus results and save them into the spike data structure
    if spikeData(ii).spikeCount>0
        clustDat = load(sprintf('%s/times_ch_%s.mat',char(folderName),char(spikeData(ii).channelName))); % Load the Wave_clus clustering data for the current channel
        if ~isempty(clustDat.cluster_class)
            spikeData(ii).units = clustDat.cluster_class(:,1); % Save the clustering data as the unit numbers for the current channel
        else
            spikeData(ii).units = zeros(size(spikeData(ii).timeStamps)); % If no units were found, save zeros instead
        end
    else
        spikeData(ii).units = []; % If no spikes were in this file, save an empty placeholder for units
    end


end

%% Save Final copy
strs = regexp(fileName,'.mcd','Split'); % Get filename from the last MCD file

finalData = []; % initialize data structure for all data
waveData = []; % initialize data structure for waveform data
for ii = 1:length(spikeData)
   timeStamps = getfield(spikeData(ii),'timeStamps'); % extract timestamp information for current channel
   units = getfield(spikeData(ii),'units'); % extract unit information for the current channel
   fCopy = [units timeStamps]; % Place the units and timestamps into single matrix
   finalData = setfield(finalData,sprintf('ch_%s',char(spikeData(ii).channelName)),fCopy); % Save units and timestamps into the finalData structure
   uniqueUnits = unique(units); % Find all the unique units in the current channel
   for jj = 1:length(uniqueUnits) % For each unique unit, save the average unit waveform into a separate data structure
       waveData = setfield(waveData,sprintf('ch_%s_unit%d',char(spikeData(ii).channelName),uniqueUnits(jj)),nanmean(spikeData(ii).waveforms(:,units==uniqueUnits(jj))')); % calculate and store average unit waveform shape
   end
   
end
finalData.waveformData = waveData; % Store waveform shape into finalData structure
finalData.A1 = A1; % Save analog 1 input into finalData structure
finalData.A2 = A2; % Save analog 2 input into finalData structure
finalData.D1 = D1; % save digital input 1 into finalData structure
finalData.D2 = D2; % save digital input 2 into finalData structure
finalData.D3 = D3; % save digital input 3 into finalData structure
finalData.maxTime = maxTime; % Save max time for the MCD recordings

% save(sprintf('%s/procFile.mat',folderName),'-struct','finalData') % Save
% a backup copy of the final data in the temproary file directory
save(saveFile,'-struct','finalData') % Save the final data
cd(oDir); % Navigate to original directory
end



function [spikeData,initialTime,dat] = readSpikeData(fileName)

%% Check for 32 or 64 bit operating system
% check = computer('arch');
% if strcmp(check,'win32')
%     DLLfileName = loader('D:\Troy\Wave_Clus UIC\nsMCDLibrary_3.7b\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary.dll',1); % Load DLL for MCD files for 32 bit systems
% else
%     DLLfileName = loader('D:\Troy\Wave_Clus UIC\nsMCDLibrary_3.7b\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary64.dll',1);% Load DLL for MCD files for 64 bit systems
% end

DLLfileName = 'D:\Troy\Wave_Clus UIC\nsMCDLibrary_3.7b\Matlab\Matlab-Import-Filter\Matlab_Interface\nsMCDLibrary64.dll';

%% Read MCD File
[nsresult] = ns_SetLibrary(DLLfileName); % Call neuroshare library to open MCD file
disp(DLLfileName)
[nsresult,info] = ns_GetLibraryInfo();
[nsresult, hfile] = ns_OpenFile(fileName); % Open current MCD file
[nsresult,info]=ns_GetFileInfo(hfile); % Get file info
for ii = 1:info.EntityCount % Go through file and extract all the data field names
   [nsresult,entity] = ns_GetEntityInfo(hfile,ii);
   labels(ii) = cellstr(entity.EntityLabel(end-1:end)); % extract data field names
end

%% Get Starting Time
refTime = datenum(2012,01,01,00,00,00); % Set reference time as January 1st 2012
timeVec = [info.Time_Year info.Time_Month info.Time_Day info.Time_Hour info.Time_Min info.Time_Sec+info.Time_MilliSec/1000]; % Get time information from MCD file
duration = info.TimeSpan; % Get timespan of MCD file
initialTimeVec = datenum(info.Time_Year,info.Time_Month,info.Time_Day,info.Time_Hour,info.Time_Min,info.Time_Sec+info.Time_MilliSec/1000); % Convert time information into a date vector
initialTime = etime(datevec(initialTimeVec),datevec(refTime)); % Calculate elapsed time between the reference and the MCD time information

%% Get all entity types from data field names
for ii = 1:length(labels) % Go through all data fields
    [nsresult,entity] = ns_GetEntityInfo(hfile,ii);
    enLab(ii) = cellstr(entity.EntityLabel); %Save whole data field name
    enType(ii) = entity.EntityType; % Save what type of entity the data field is
    itemCount(ii) = entity.ItemCount; % Save how many items are in the entity
end


%% Get Spike Data

spikeInd = find(strncmp(enLab,'spks',4)); % Find all data entities that correspond to spikes
spikeLab = enLab(spikeInd);
spikeType = enType(spikeInd);
spikeCount = itemCount(spikeInd);

for ii = 1:length(spikeInd)
    channelStr(ii) = cellstr(spikeLab{ii}(end-1:end)); % Extract the channel name for the current spike data
    channelNum(ii) = str2num(char(channelStr(ii))); % Extract the channel number for the current spike data
end

%% Get all spike info

for ii = 1:length(spikeInd) % Go through MCD file and extract and save spike data to data structure
    entityID = spikeInd(ii); % get current entity ID
    spikeData(ii).channelName = channelStr(ii); % save channel name
    spikeData(ii).channelNumber = channelNum(ii); % save channel number
    spikeData(ii).spikeCount = spikeCount(ii); % save the number of spikes
    timeStamps = []; % initialize the timestamps array
    waveforms = []; % intialize the waveforms array
    if spikeCount(ii) > 0 % If there are spikes in the current channel, extract them and their timestamps
        [ns_RESULT, timeStamps, waveforms, SampleCount, UnitID] = ns_GetSegmentData(hfile, entityID, 1:spikeCount(ii));
    end
    spikeData(ii).timeStamps = timeStamps+initialTime; % Save timestamps with actual real-world times
    spikeData(ii).waveforms = waveforms; % Save spike waveform shapes
end


%% Get Digital Data
entitityID = find(strcmp(labels,'D1')==1,1,'first'); % Find the entity for digital data
if ~isempty(entitityID) % If the MCD file has digital data
    [nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
    [nsresult,a] = ns_GetAnalogInfo(hfile,entitityID); 
    [nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount); % Extract all digital data
    d1 = bitget(data,1); % Save digital data as binary array for digital input 1
    d2 = bitget(data,2); % Save digital data as binary array for digital input 2
    d3 = bitget(data,3); % Save digital data as binary array for digital input 3
else
    d1 = [];
    d2 = [];
    d3 = [];
end

%% Get Analog Data
entitityID = find(strcmp(labels,'A1')==1,1,'first'); % Find the entity for analog input 1
if ~isempty(entitityID)
    [nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
    [nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
    [nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount); % Extract data
else
    data = [];
end
a1 = data; % Save analog input 1

entitityID = find(strcmp(labels,'A2')==1,1,'first');
if ~isempty(entitityID)
[nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
[nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
[nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount); % Extract data
else
    data = [];
end
a2 = data; % Save analog input 2

entitityID = find(strcmp(labels,'A3')==1,1,'first');
if ~isempty(entitityID)
[nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
[nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
[nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount); % Extract data
else
    data = [];
end
a3 = data; % Save analog input 3


%% Get time

[nsresult,info]=ns_GetFileInfo(hfile);
% dat.finalTime = info.TimeSpan+initialTime;
finalTimeVec = addtodate(initialTimeVec,info.TimeSpan*1000,'millisecond'); % Get the time when the MCD file ended in milliseconds
dat.finalTime = etime(datevec(finalTimeVec),datevec(refTime)); % Get real-world time for the final time
dat.initialTime = initialTime; % Save starting time for MCD file
time = (info.TimeStampResolution:info.TimeStampResolution:info.TimeSpan)';


%% Get Stimulus, Cycle, and Pulse Tags


dat.D1 = time(find(diff(d1)>0)+1)+initialTime; % Find and save all rising edges from digital input 1
dat.D2 = time(find(diff(d2)>0)+1)+initialTime; % Find and save all rising edges from digital input 2
dat.D3 = time(find(diff(d3)>0)+1)+initialTime; % Find and save all rising edges from digital input 3

%% Get rising and falling edges of a1 and a2
% Light Waveforms
try
da1 = diff(a1)./diff(time); % Calculate the derivative of the analog input 1
[pks1,risingEdges] = findpeaks(da1,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time)))); % Find rising edges of data
risingEdges = time(risingEdges+1); % Get time points for rising edges
[pks2,fallingEdges] = findpeaks(-da1,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time)))); % Find falling edges of data
fallingEdges = time(fallingEdges+1); % Get time points for falling edges
dat.A1 = [risingEdges+initialTime; -(fallingEdges+initialTime)]; % Save rising and falling edges into data structure as positive and negative times, respectively
catch
    dat.A1 = [];
end

% Injection Edges
try
da2 = diff(a2)./diff(time);% Calculate the derivative of the analog input 2
[pks1,risingEdges] = findpeaks(da2,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time)))); % Find rising edges of data
risingEdges = time(risingEdges+1); % Get time points for rising edges
[pks2,fallingEdges] = findpeaks(-da2,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time)))); % Find falling edges of data
fallingEdges = time(fallingEdges+1); % Get time points for falling edges
dat.A2 = [risingEdges+initialTime; -(fallingEdges+initialTime)]; % Save rising and falling edges into data structure as positive and negative times, respectively
clear time
catch
    dat.A2 = [];

end





end

function newData = appendSpikeData(spikeData,cSpikeData)
% Subfunction for saving spike data into data structure

if isempty(spikeData) % If this is the first time the subfunction is used
    for ii = 1:length(cSpikeData)        
        [m,n]=size(cSpikeData(ii).waveforms); % Get number of spikes
        test = find([m n]==30); % find the orientation of the spike data
        newData(ii).spikeCount = length(cSpikeData(ii).timeStamps); % count the spikes
        newData(ii).timeStamps = [cSpikeData(ii).timeStamps]; % save all spike timestamps
        if test==2
            newData(ii).waveforms = [cSpikeData(ii).waveforms']; % save spike waveforms in the appropriate orientation
        else
            newData(ii).waveforms = [cSpikeData(ii).waveforms]; % save spike waveforms in the appropriate orientation
        end
        newData(ii).channelName = cSpikeData(ii).channelName; % save channel name
        newData(ii).channelNumber = cSpikeData(ii).channelNumber; % save channel number
    end
else % for subsequent calls to this function
    newChannels = cell2mat(arrayfun(@(x) cSpikeData(x).channelNumber,1:length(cSpikeData),'Unif',false))'; % get all channel names from the new MCD file
    for ii = 1:length(spikeData)
        if ~isempty(spikeData(ii).channelNumber)
            ind = find(spikeData(ii).channelNumber==newChannels); % Find the index of the current channel in the newChannels array
        else
            ind = [];
        end
%         if isempty(ind)
%             error('Missing Data!')
%         end
        if ~isempty(ind)
            newData(ii).spikeCount = spikeData(ii).spikeCount+length(cSpikeData(ind).timeStamps); % Append spike count from new MCD file onto same channel as previous files
            newData(ii).timeStamps = [spikeData(ii).timeStamps; cSpikeData(ind).timeStamps]; % Append time stamps
            [m,n]=size(cSpikeData(ind).waveforms); % Get size of current spike waveforms
            test = find([m n]==30);
            if test==2
                newData(ii).waveforms = [spikeData(ii).waveforms cSpikeData(ind).waveforms']; % Append current waveforms
            else
                newData(ii).waveforms = [spikeData(ii).waveforms cSpikeData(ind).waveforms]; % Append current waveforms
            end
            newData(ii).channelName = spikeData(ii).channelName; % Save channel Name
            newData(ii).channelNumber = spikeData(ii).channelNumber; % save channel number
        end
    end
    
end





end
