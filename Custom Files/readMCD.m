function varargout = readMCD(fileName,flag,spikeFlag,rawElectrodeFlag,a3Flag,tFlag)
%Reads a provided MCD file and returns the rising (positive) and falling
%(negative) edges of light (A1) and injection (A2) outputs. Also outputs
%the stimulus tags (D1), cycle tags (D2), pulse tags (D3), and time
%extracted from the MCD file
rawElectrodeData=[];
%% Check for 32 or 64 bit operating system
check = computer('arch');
if strcmp(check,'win32')
    DLLfileName = loader('nsMCDlibrary.dll',1);
else
    DLLfileName = loader('nsMCDlibrary64.dll',1);
end

%% Check for flags

if nargin == 1
    flag = 0;
    spikeFlag = false;
    rawElectrodeFlag = false;
    a3Flag = false;
    tFlag = false;
end

if nargin == 2
    spikeFlag = false;
    rawElectrodeFlag = false;
    a3Flag = false;
    tFlag = false;
end

if nargin == 3
    rawElectrodeFlag = false;
    a3Flag = false;
    tFlag = false;
end

if nargin == 4
    a3Flag = false;
    tFlag = false;
end

if nargin == 5
    tFlag = false;
end

%% Read File
[nsresult] = ns_SetLibrary(DLLfileName);
[nsresult,info] = ns_GetLibraryInfo();
[nsresult, hfile] = ns_OpenFile(fileName);
[nsresult,info]=ns_GetFileInfo(hfile);
for ii = 1:info.EntityCount
   [nsresult,entity] = ns_GetEntityInfo(hfile,ii);
   labels(ii) = cellstr(entity.EntityLabel(end-1:end)); 
   types(ii) = cellstr(entity.EntityLabel(1:4));
end

%% Get Digital Data
entitityID = find(strcmp(labels,'D1')==1,1,'first');
if ~isempty(entitityID)
    [nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
    [nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
    [nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount);
    d1 = bitget(data,1);
    d2 = bitget(data,2);
    d3 = bitget(data,3);
else
    d1 = [];
    d2 = [];
    d3 = [];
end



%% Get Analog Data
entitityID = find(strcmp(labels,'A1')==1,1,'first');
if ~isempty(entitityID)
    [nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
    [nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
    [nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount);
else
    data = [];
end
a1 = data;

entitityID = find(strcmp(labels,'A2')==1,1,'first');
if ~isempty(entitityID)
[nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
[nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
[nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount);
else
    data = [];
end
a2 = data;

entitityID = find(strcmp(labels,'A3')==1,1,'first');
if ~isempty(entitityID)
[nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
[nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
[nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount);
else
    data = [];
end
a3 = data;


%% Get Spike Data

if spikeFlag
    spikeInd = find(strcmp(types,'spks'));  
    for ii = 1:length(spikeInd)
        entityID = spikeInd(ii);
        [nsresult,entity] = ns_GetEntityInfo(hfile,entityID);
        spikeCount = entity.ItemCount;
        if spikeCount > 0
            [ns_RESULT, timeStamps, waveforms, SampleCount, UnitID] = ns_GetSegmentData(hfile, entityID, 1:spikeCount);
            timeStamps = round(timeStamps*10000)/10000;
            waveforms = waveforms';
        else
            timeStamps = [];
            waveforms = [];
        end
        spikeData(ii).channel = str2num(char(labels(entityID)));
        spikeData(ii).channelName = strcat('ch_',char(labels(entityID)));
        spikeData(ii).spikeCount = spikeCount;
        spikeData(ii).timeStamps = timeStamps;
        spikeData(ii).waveforms = waveforms;
        spikeData(ii).entityInfo = entity;
    end
end


%% Get Raw Electrode Data

if rawElectrodeFlag
    elecInd = find(strcmp(types,'elec'));  
    for ii = 1:length(elecInd)
        entityID = elecInd(ii);
        [nsresult,entity] = ns_GetEntityInfo(hfile,entityID);
        sampleCount = entity.ItemCount;
        [nsresult,a] = ns_GetAnalogInfo(hfile,entityID);
        [nsresult,count,data]=ns_GetAnalogData(hfile,entityID,1,entity.ItemCount);
        rawElectrodeData(ii).channel = str2num(char(labels(entityID)));
        rawElectrodeData(ii).channelName = strcat('ch_',char(labels(entityID)));
        rawElectrodeData(ii).sampleCount = sampleCount;
        rawElectrodeData(ii).data = data;
        rawElectrodeData(ii).analogInfo = a;
        rawElectrodeData(ii).entityInfo = entity;
    end
else 
    rawElectrodeData = [];
end
    

%% Get time
try
entitityID = find(strcmp(labels,'A1')==1,1,'first');
[nsresult,entity] = ns_GetEntityInfo(hfile,entitityID);
[nsresult,a] = ns_GetAnalogInfo(hfile,entitityID);
[nsresult,count,data]=ns_GetAnalogData(hfile,entitityID,1,entity.ItemCount);
catch ME
    fprintf('%s',ME.message)
end

try
    [ns_RESULT, time] = ns_GetTimeByIndex(hfile, 1, 1:length(data));
    if length(time)~=length(data)
        error('Something wrong with time vector')
    end
catch
    clear time
    try
    for ii = 1:count
        [ns_RESULT, time(ii)] = ns_GetTimeByIndex(hfile, 1, ii);
    end
    catch
        clear time
        [nsresult,info]=ns_GetFileInfo(hfile);
        time = (info.TimeStampResolution:info.TimeStampResolution:info.TimeSpan)';
        if length(time) ~= count
            error('Wrong time stamp values')
        end
    end
end

%% Get Starting Time
refTime = datenum(2012,01,01,00,00,00);
timeVec = [info.Time_Year info.Time_Month info.Time_Day info.Time_Hour info.Time_Min info.Time_Sec+info.Time_MilliSec/1000];
duration = info.TimeSpan;
initialTimeVec = datenum(info.Time_Year,info.Time_Month,info.Time_Day,info.Time_Hour,info.Time_Min,info.Time_Sec+info.Time_MilliSec/1000);
initialTime = etime(datevec(initialTimeVec),datevec(refTime));

%% Get Stimulus, Cycle, and Pulse Tags
D1 = time(find(diff(d1)>0)+1);
D2 = time(find(diff(d2)>0)+1);
D3 = time(find(diff(d3)>0)+1);

%% Get rising and falling edges of a1 and a2
% Light Waveforms
da1 = diff(a1)./diff(time);
[pks1,risingEdges] = findpeaks(da1,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time))));
risingEdges = time(risingEdges+1);
[pks2,fallingEdges] = findpeaks(-da1,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time))));
fallingEdges = time(fallingEdges+1);
A1 = [risingEdges; -fallingEdges];

% Injection Edges
da2 = diff(a2)./diff(time);
[pks1,risingEdges] = findpeaks(da2,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time))));
risingEdges = time(risingEdges+1);
[pks2,fallingEdges] = findpeaks(-da2,'minpeakheight',500,'minpeakdistance',round(0.001/mean(diff(time))));
fallingEdges = time(fallingEdges+1);
A2 = [risingEdges; -fallingEdges];



% %% Get Individual Cycle Waveforms
% 
% for ii = 1:length(D1)
%     t1 = D1(ii);
%     if ii == length(D1)
%         t2 = time(end);
%     else
%         t2 = D1(ii+1);
%     end
%     
%     pulseTags = D3(find(D3>=t1&D3<t2));
% 
%     cycleTags = D2(find(D2>=t1&D2<t2));
%     cycleLength = mean(diff(cycleTags));
%     limitLength = round(cycleLength/mean(diff(time)))-20;
%     if isnan(limitLength)
%         limitLength = 0;
%     end
%     lightWFs = nan(length(cycleTags),limitLength);
%     injectionWFs = nan(length(cycleTags),limitLength);    
%     for jj = 1:length(cycleTags)
%        c1 = cycleTags(jj);
%        if jj == length(cycleTags)
%            c2 = c1+cycleLength;
%        else
%            c2 = cycleTags(jj+1);
%        end
%        indices = find(time >= c1 & time < c2);
%        if length(indices) >= limitLength
%             lightWFs(jj,:) = a1(indices(1:limitLength));
%             injectionWFs(jj,:) = a2(indices(1:limitLength));
%        else
%             lightWFs(jj,1:length(indices)) = a1(indices);
%             injectionWFs(jj,1:length(indices)) = a2(indices);
%        end
%         
%     end
%     lightWF = nanmean(lightWFs);
%     injectionWF = nanmean(injectionWFs);
%     timeWF = linspace(0,cycleLength,limitLength);
%     nt = 1;
% 
% 
% end
[nsresult,info]=ns_GetFileInfo(hfile);
% dat.finalTime = info.TimeSpan+initialTime;
finalTimeVec = addtodate(initialTimeVec,info.TimeSpan*1000,'millisecond');
finalTime = etime(datevec(finalTimeVec),datevec(refTime));

if flag == 1
    A1 = a1;
    A2 = a2;
end
if nargout == 6
    varargout{1} = A1;
    varargout{2} = A2;
    varargout{3} = D1;
    varargout{4} = D2;
    varargout{5} = D3;
    varargout{6} = time;
end

if nargout == 7
    varargout{1} = A1;
    varargout{2} = A2;
    varargout{3} = D1;
    varargout{4} = D2;
    varargout{5} = D3;
    varargout{6} = time;
    
    if spikeFlag
        varargout{7} = spikeData;
    elseif flag
        varargout{7} = initialTime;
    elseif a3Flag
        varargout{7} = a3;
    else
        varargout{7} = rawElectrodeData;
    end

    
end


if nargout == 8
    varargout{1} = A1;
    varargout{2} = A2;
    varargout{3} = D1;
    varargout{4} = D2;
    varargout{5} = D3;
    varargout{6} = time;
    varargout{8} = rawElectrodeData;
    
    if spikeFlag
        varargout{7} = spikeData;
    elseif flag
        varargout{7} = initialTime;
        varargout{8} = finalTime;
    elseif a3Flag
        varargout{7} = a3;
    elseif tFlag
        varargout{7} =[initialTime finalTime];
    else
        varargout{7} = spikeData;
    end
  
end











