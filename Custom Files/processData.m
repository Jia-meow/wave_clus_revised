function processData(varargin)
% Processes data and saves the requisite information to the master mat file
% and other places
%% Get Inputs
fileName = varargin{1};
if nargin == 2
    masterFlag = varargin{2};
else
    masterFlag = [];
end
flagger = 0;

%% Append file if manually sorted

if ischar(masterFlag)|iscell(masterFlag)
   originalFile = load(char(masterFlag));
   newFile = load(char(fileName));
%    originalFields = fieldnames(originalFile);
%    newFields = fieldnames(newFile);
   originalFields = fields(originalFile);
   newFields = fields(newFile);
   
   ind = find(~ismember(originalFields,newFields)&~strcmp(originalFields,'time'));
   
   for ii = 1:length(ind)
       newFile = setfield(newFile,char(originalFields(ind(ii))),getfield(originalFile,char(originalFields(ind(ii)))));
   end
   newFile.maxTime = max(originalFile.time);
   newFile.minTime = min(originalFile.time);
   save(char(fileName),'-struct','newFile');
   masterFlag = []; 
end

%% Check for batch process

if ischar(fileName)
    
%% Get date and recording name
startInd = findstr(fileName,'MEA')+4;
if length(startInd)>1
    startInd = startInd(end);
end
dateInd = find(fileName(startInd:end)=='/',1,'first')-2+startInd;
recordingInd = find(fileName(startInd:end)=='/',1,'last')-2+startInd;
date = fileName(startInd:dateInd);
recordingName = fileName(dateInd+2:recordingInd);
fprintf('%s - %s\n',date,recordingName)
fprintf('Starting Processing...\n');


%% Get raw data
rawData = load(fileName);

%% Check to see if MCD data has been incorporated and load if not
if ~isfield(rawData,'D1')
    preProcessData(fileName);
    rawData = load(fileName);
end


%% Determine which substructures contain channel info
% fieldNames = fieldnames(rawData);
% channelNames = fieldNames(strncmp(fieldNames,'ch_',3));
fieldNames = fields(rawData);
channelNames = fieldNames(strncmp(fieldNames,'ch_',3));
if isempty(channelNames)
    error('No Channels Listed in MAT File!')
end

%% Set number of spikes threshold
frequencyThreshold = 0.3; % frequency in Hz ... only process a channel if channel spike frequency is greater than this number
if sum(strcmp(fieldNames,'maxTime'))==0
    maxRecordingTime = rawData.time(end);
else
    maxRecordingTime = rawData.maxTime;
end
numberOfSpikesThreshold = frequencyThreshold * maxRecordingTime;

%% Get channels of interest
for ii = 1:length(channelNames)
    channelInfo = getfield(rawData,char(channelNames(ii)));
    if ~isempty(channelInfo)
        waveformNumbers(ii) = length(unique(channelInfo(:,2)));
        unitNumbers(ii) = sum(unique(channelInfo(:,1))~=0);
    else
        waveformNumbers(ii) = 0;
        unitNumbers(ii) = 0;
    end
end
goodChannelsInd = find(waveformNumbers>numberOfSpikesThreshold); %find channels which have firing frequency > frequencyThreshold
channelsOfInterest = channelNames(goodChannelsInd);

%% Get Stimulus Information from Stimulus File
XMLfile = strcat(fileName(1:find(fileName=='/',1,'last')),'stimulusFile.xml');
if ~exist(XMLfile,'file')
    restructureStimulusFile(fileName(1:find(fileName=='/',1,'last')),0);
    flagger = 1;
end
stimulusFile = readStimulusFile(XMLfile);
stimNames = arrayfun(@(x) stimulusFile(x).name,1:length(stimulusFile),'Unif',false)';
if sum(strncmp(stimNames,'Spot S',6))>0 % Checks to see if there are any spot stimulus visual stimuli
    D1 = rawData.D1;
    D3 = rawData.D3;
    newD1 = [];
    d1Counter = 1;
    lastFlag = 0;
    for ii = 1:length(stimulusFile)
        if ~strncmp(stimulusFile(ii).name,'Spot S',6)
            if lastFlag == 3
                newD1(ii) = newD1(end)+str2num(stimulusFile(ii-1).controlVals{1})*stimulusFile(ii-1).cycleNum+2;
                lastFlag = 0;
            else
                newD1(ii) = D1(d1Counter);
                d1Counter = d1Counter+1;
                lastFlag = 1;
            end
        else
            i1 = newD1(end);
            d3i = find(D3>i1,1,'first');
            newD1(ii) = D3(d3i);
            lastFlag = 3;
        end
    end
    newD1(end+1) = D1(end);
    D1 = newD1;
    D2 = rawData.D2;
    newD2 = [];
    for ii = 1:length(stimulusFile)
        s1 = D1(ii);
        s2 = D1(ii+1);
        if strncmp(stimulusFile(ii).name,'Spot Sp',7)
            cD2 = s1+0.1:str2num(stimulusFile(ii).controlVals{1}):s2;
            cD2 = cD2(1:stimulusFile(ii).cycleNum)';
        else            
            cD2 = D2(D2>=s1 &D2<=s2);
        end
        newD2 = [newD2; cD2];
    end
            
    rawData.D2 = newD2;
    rawData.D1 = newD1;
end

numberOfStimuli = length(rawData.D1)-1;
if length(stimulusFile)< numberOfStimuli & flagger == 0
    restructureStimulusFile(fileName(1:find(fileName=='/',1,'last')),1);
end
if length(stimulusFile)< numberOfStimuli & flagger == 1
    restructureStimulusFile(fileName(1:find(fileName=='/',1,'last')),2);
    stimulusFile = readStimulusFile(XMLfile);
    if length(stimulusFile)< numberOfStimuli
        error('Missing stimuli!')
    end
end
%% Query Master for location of where to place processed data
masterFileName = loader('MEA/master.mat',1);
master = load(masterFileName);
if isempty(master.processedDataFolder)
    processedDataFolder = strcat(fileName(1:find(fileName=='/',1,'last')),'processedData');
else
    startInd = strfind(fileName,'MEA')+3;
    if length(startInd)>1
        startInd = startInd(end);
    end
    endInd = find(fileName == '/',1,'last')-1;
    processedDataFolder = strcat(master.processedDataFolder,fileName(startInd:endInd));
end

%% Create folder for processed data if one does not exist
dirExistence = dir(processedDataFolder);
if isempty(dirExistence)
   mkdir(processedDataFolder); 
end

%% Process and save the data
numberOfStimuli = length(rawData.D1)-1;
N = sum(unitNumbers(goodChannelsInd));
reverseStr = '';
unitNumber = 1;
for ii = 1:length(channelsOfInterest)
    channelInfo = getfield(rawData,char(channelsOfInterest(ii)));
    currentUnits = unique(channelInfo(:,1));
    currentUnits = currentUnits(currentUnits~=0);
    
    for jj = 1:length(currentUnits); % Go through units of this channel and then extract/save information to separate .mat file
        unitInd = channelInfo(:,1) == currentUnits(jj);
        allSpikes = channelInfo(unitInd,2);
        clear unit
        unit.channelName = char(channelsOfInterest(ii));
        unit.channelNumber = str2num(unit.channelName(end-1:end));
        unit.unitNumber = unitNumber;
        unit.channelUnitNumber = currentUnits(jj);
        for kk = 1:numberOfStimuli
            unit.response(kk).stimulusName = stimulusFile(kk).name;
            s1 = rawData.D1(kk);
            s2 = rawData.D1(kk+1);
            spikes = allSpikes(allSpikes >= s1 & allSpikes <= s2);
            cycleTags = rawData.D2(rawData.D2>=s1 &rawData.D2<=s2);
            try
                [kern, kernTime,rasterX,rasterY] = kernelZero(spikes,cycleTags); % Try to calculate gaussian kernel density estimation of spike rate
            catch
                kern = NaN;
                kernTime = NaN;
                rasterX = NaN;
                rasterY = NaN;
            end
            unit.response(kk).kde = kern;
            unit.response(kk).kdeTime = kernTime;
            unit.response(kk).rasterX = rasterX;
            unit.response(kk).rasterY = rasterY;
            unit.response(kk).cycleTags = cycleTags;
            
            
        end
        processedFileName = sprintf('%s/unit_%0.3d',processedDataFolder,unitNumber);
        save(processedFileName,'unit') % Save data for current unit
        msg = sprintf('Processed %d/%d', unitNumber, N);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        unitNumber = unitNumber + 1;

    end

end
fprintf('\n')
%% Update Master File

try
    if isempty(masterFlag)
        updateMaster(fileName)
    else
        updateMaster(fileName,masterFlag);
    end
    fprintf('Master Updated...\n')
catch ME
    if strcmp(ME.message,'Redundant')
       fprintf('Already in Master file...\n') 
    end
end


%% Save spike shapes
try
    extractUnitSpikeShapes(sprintf('%s/%s',date,recordingName),false); % Old file that doesn't work anymore- ignore
catch
    folderName = sprintf('%s/%s',date,recordingName);
    str = regexp(folderName,'ing','split');
    recNum = str2num(str{2});
    matFile = loader(sprintf('%s/r%d.mat',folderName,recNum));
    matDat = load(matFile);
    [unitFiles,units] = getUnitFiles(folderName);
    for ii = 1:length(unitFiles)
        cUnit = load(char(unitFiles(ii)));
        chList(ii) = cellstr(cUnit.unit.channelName);
        chNum(ii) = cUnit.unit.channelNumber;
        chUnitNum(ii) = cUnit.unit.channelUnitNumber;
        waveData = getfield(matDat.waveformData,char(sprintf('%s_unit%d',char(chList(ii)),chUnitNum(ii))));
        cUnit = setfield(cUnit.unit,'spikeShape',waveData);
        unit = cUnit;
        save(char(unitFiles(ii)),'unit')
    end
end
    
    
%% Batch Job
else
   clc
   for pp = 1:length(fileName)
       currentFile = char(fileName(pp));
       startInd = findstr(currentFile,'MEA')+4;
       if length(startInd)>1
           startInd = startInd(end);
       end
       dateInd = find(currentFile(startInd:end)=='/',1,'first')-2+startInd;
       recordingInd = find(currentFile(startInd:end)=='/',1,'last')-2+startInd;
       date = currentFile(startInd:dateInd);
       recordingName = currentFile(dateInd+2:recordingInd);
%        try
           fprintf('File %d/%d: ',pp,length(fileName))
            processData(currentFile,masterFlag);
%        catch ME
%             fprintf('Error: %s\n',ME.message)
%        end
          
   end
    
end

end

function [varargout] = kernelZero(varargin)


spikes = varargin{1};
cycleTags = varargin{2};

    compr = ceil(length(min(cycleTags):0.0001:max(spikes))/4000); % Finds compression factor for gaussian KDE
    numberOfCycles = length(cycleTags); % Find number of cycles
    rasterX = [];
    rasterY = [];
        cycleTimer = mean(diff(cycleTags)); % Get the trial/cycle time
        for jj = 1:numberOfCycles
                t1 = cycleTags(jj); % Get start time of current trial
                t2 = (cycleTags(jj)+cycleTimer); % Get end time of current trial
                cSpikes = spikes(find(spikes>=t1 & spikes<=t2)); % Extract spikes of current trial
                rasterX = [rasterX (cSpikes'-cycleTags(jj))]; % Save timestamps of current trial in a format for rasters and KDE
                rasterY = [rasterY (ones(size(cSpikes'))*jj)]; % save trial information for current timestamps
        end
        [kern,kernTime] = gausskernZero(75,spikes,min(cycleTags),max(spikes),0.0001,compr); % Calculate gaussian KDE with a sigma of 75
        for hh = 1:length(cycleTags)
            t1 = cycleTags(hh);
            t2 = cycleTags(hh)+cycleTimer;
            ind = find(kernTime >= t1 & kernTime <= t2); % find KDE for given trial
            cycleKern{hh} = kern(ind); % save trial KDE into array
            cycleTime{hh} = kernTime(ind)-cycleTags(hh); % Save trial times into array
            len(hh) = length(kernTime(ind));
        end
        ind = len==0;
        len(ind) = NaN;
        maxlen = floor(nanmedian(len)); % Get median length of trial KDEs since they may be different lengths
        avgKern = nan(length(cycleTime),maxlen); % create data array to hold trial KDEs
        avgKernTime = nan(length(cycleTime),maxlen);
        for hh = 1:length(cycleTime)
            if length(cycleKern{hh})>= maxlen % If current trial KDE is longer than the median length, it will be truncated
                avgKern(hh,1:maxlen) = cycleKern{hh}(1:maxlen);
                avgKernTime(hh,1:maxlen) = cycleTime{hh}(1:maxlen);
            elseif ~isempty(cycleKern{hh}) %Otherwise it will all be saved. if it is shorter than the median length it is padded with NaNs
                avgKern(hh,1:length(cycleKern{hh})) = cycleKern{hh};
                avgKernTime(hh,1:length(cycleKern{hh})) = cycleTime{hh};       
            end
        end
        avgKern = nanmean(avgKern); % Get average KDE
        avgKernTime = nanmean(avgKernTime); % Get average timing for KDE
        avgKernTime = linspace(min(avgKernTime),max(avgKernTime),length(avgKern)); % Ensure the average timing is consistent
        varargout{1} = avgKern;
        varargout{2} = avgKernTime;
        if nargout == 4
            varargout{3} = rasterX;
            varargout{4} = rasterY;
        end
    
        



end

function [final_signal,time] = gausskernZero(sigma,timestamps,start,stop,dt, compr)


% dt = floor(min(diff(timestamps))*1000)/1000;
% timestamps = timestamps - start;
bin_edges = start:dt:stop;
spiketrain = histc(timestamps,bin_edges); % Convert spike timestamps into peristimulus histogram

gauss_duration = 3*sigma*2; % Set up gaussian window
gauss_t = -gauss_duration/2:dt*1000:gauss_duration/2;
gauss = 1/sqrt(1*pi*sigma^2)*exp(-gauss_t.^2/(2*sigma^2));

final_signal = conv(double(spiketrain),gauss,'same')*1000; % Perform convolution with guassian window
time = bin_edges;
ind = 1:compr:length(final_signal); % Take every "compr" samples to reduce size of resulting gaussian KDE
final_signal = final_signal(ind);
time = time(ind);
end