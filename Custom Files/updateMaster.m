function updateMaster(varargin)
% Function will open the master mat file and attempt to update it with the
% date and recording information contained in fileName (varargin{1}). If
% you want to force an update or alter some other part of the function
% include the following flags as varargin{2}:
% flag = 1 --> Force update of master.mat. Will overwrite any information
% currently contained in master.mat for the given date and recording

%% Input variables

fileName = varargin{1};
if nargin >1
    flag = varargin{2};
else
    flag = 0;
end

%% Initialization
masterFileName = loader('MEA/master.mat',1);
try
    master = load(masterFileName);
catch
    error('Master not present!');
end

%% Check if date and recording are listed and exit if they need be
startInd = findstr(fileName,'MEA')+4;
if length(startInd)>1
    startInd = startInd(end);
end
dateInd = find(fileName(startInd:end)=='/',1,'first')-2+startInd;
recordingInd = find(fileName(startInd:end)=='/',1,'last')-2+startInd;
date = fileName(startInd:dateInd);
recordingName = fileName(dateInd+2:recordingInd);
recordingNumber = str2double(recordingName(find(recordingName == 'g',1,'last')+1:end));

if sum(strcmp(master.date,date)&strcmp(master.recordingName,recordingName))>0&flag~=1
    error('Redundant')
end

%% update information in master.mat

allDates = master.date;
allRecordingNames = master.recordingName;
allRecordingNumbers = master.recordingNumber;

if flag == 1
    ii = find(strcmp(master.date,date)&strcmp(master.recordingName,recordingName)==1,1,'first');
else
    ii = length(allDates)+1;
end

allDates(ii) = cellstr(date);
allRecordingNames(ii) = cellstr(recordingName);
allRecordingNumbers(ii) = recordingNumber;


%% Update Stimulus Info
stimuli = master.stimulusInfo;
stimulusFileName = loader(sprintf('MEA/%s/%s/stimulusFile.xml',date,recordingName),1);
stimuli{ii} = readStimulusFile(stimulusFileName);


%% Resort the dates, recordings, etc

sortInd = reSort(allDates,allRecordingNumbers);
master.date = allDates(sortInd);
master.recordingName = allRecordingNames(sortInd);
master.recordingNumber = allRecordingNumbers(sortInd);
master.stimulusInfo = stimuli(sortInd);

master.allParameters = gatherMasterData(master);

%% Save master file
save(masterFileName,'-struct','master')

end

function sortInd = reSort(dates,recordingNums)


[r1,i1] = sort(recordingNums);
d1 = dates(i1);
[d2,i2] = sort(d1);
sortInd = (i1(i2));

end

function newDat = gatherMasterData(master)
% Resort the master data and extract relevant information from each one for
% indexing purposes

[stimulusTypes,stimTypeNames] = readStimulusTypes();

counter =1;
for ii = 1:length(master.stimulusInfo)
    
    for jj = 1:length(master.stimulusInfo{ii})
        newDat.date(counter) = str2num(char(master.date{ii}));
        newDat.recording(counter) = cellstr(master.recordingName(ii));
        newDat.name(counter) = cellstr(master.stimulusInfo{ii}(jj).name);
        newDat.allFiles(counter) = cellstr(sprintf('%s/%s',char(master.date{ii}),char(master.recordingName{ii})));
        newDat.allStims(counter) = jj;    
        newDat.chemical(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'chemical',1);
        newDat.concentration(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'concentration',0);
        newDat.pressure(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'pressure',0);
        newDat.depth(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'depth',0);
        loc = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'location',0);
        if length(loc) > 1
            [targetPixel,um,electrodes] = averageTargetPixel(sprintf('%s/%s',char(master.date{ii}),char(master.recordingName(ii))));
            [dists] = dist(targetPixel,loc);
            [m,i] = min(dists*25*um);
            newLoc = electrodes(i)+round(m)/1000;
            newDat.location(counter) = newLoc;
        else
            newDat.location(counter) = loc;
        end
        newDat.tipDiameter(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'tipDiameter',0);
        newDat.electrodeDistance(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'electrodeDistance',0);
        newDat.orientation(counter) = getParamValue(master.stimulusInfo{ii}(jj).parameterVals,master.stimulusInfo{ii}(jj).parameterNames,'orientation',1);
        ind = strcmp(stimTypeNames,newDat.name(counter));
        if sum(strcmp(stimulusTypes(ind).command,'inject'))>0
            if sum(strcmp(stimulusTypes(ind).command,'setInjectionTime'))>0
                str = 'setInjectionTime';
                
            else
                str = 'inject';
            end
            newDat.injectionTime(counter) = getParamValue(master.stimulusInfo{ii}(jj).controlVals,master.stimulusInfo{ii}(jj).controlNames,str,0);
        else
            newDat.injectionTime(counter) = NaN;
        end
        counter = counter + 1;
    end
    
end
end


function output = getParamValue(cParamVals,cParamNames,cName,type)
% Extract current parameter from the stimulus files of recordings

if type == 1
    output = cellstr(cParamVals(strcmp(cParamNames,cName)));
    if isempty(output)
        output = {''};
    end
else
    output = str2num(char(cParamVals(strcmp(cParamNames,cName))));
    if isempty(output)
        output = NaN;
    end
end
end

