function restructureStimulusFile(folderName,flag)
% Load stimulus file and make sure its in the correct format. If multiple
% stimuli were used, this combines them into a single stimulus file

if ~strcmp(folderName(end),'/')
    fN = strcat(folderName,'/');
else
    fN = folderName;
end

dirNames = dir(folderName);

for ii = 1:length(dirNames)
    cName = dirNames(ii).name;
    fileNames(ii) = cellstr(cName);
    if length(cName)>8
        if strcmp('Log',cName(end-6:end-4))
            logInd(ii) = true;
        else
            logInd(ii) = false;
        end
        if strcmp('xml',cName(end-2:end))
            xmlInd(ii) = true;
        else
            xmlInd(ii) = false;
        end
    else
        logInd(ii) = false;
        xmlInd(ii) = false;
    end
end
stimInd = strncmp(fileNames,'stim',4)&~strcmp(fileNames,'stimulusParameterInfo.mat');

%% Test for other types of stimulus files
irregularStims = fileNames(~xmlInd&~logInd&stimInd)';

%% set up buffer
buffer.name.Text = 'Empty';
buffer.cycles.Text = '1';
buffer.control.wait{1}.Text = 0.01;

%% Order Stimulus Files
% stimulusFiles = fileNames(xmlInd&~logInd&stimInd)';
stimulusFiles = fileNames(~logInd&stimInd)';

if flag == 0|flag==2
    stimulusFiles = stimulusFiles(~strcmp(stimulusFiles,'stimulusFile.xml'));
end

for ii = 1:length(stimulusFiles)
    cFile1 = isstrprop(stimulusFiles{ii},'digit');
    i1 = find(cFile1,1,'first');
    i2 = find(cFile1==0&1:length(cFile1)>i1,1,'first')-1;
    stimFileNums(ii) = str2num(char(stimulusFiles{ii}(i1:i2)));
end
[y,i] = sort(stimFileNums);
stimulusFiles = stimulusFiles(i);

%% load and make new structure


stimFile = struct([]);
for ii = 1:length(stimulusFiles)
    endPoints = 1;
    cFileName = char(stimulusFiles(ii));
    if sum(ismember(irregularStims,cFileName))>0
        text = fileread(strcat(fN,cFileName));
        ind = strfind(text,'chNum');
        i1 = ind(1);
        ind = strfind(text,'spontaneous');
        i12 = ind(1);
        ind = strfind(text,';');
        i2 = ind(find(ind>i12,1,'first'));
        eval(text(i1:i2))
        
        ind = strfind(text,'background');
        i12 = ind(1);
        ind = strfind(text,';');
        i2 = ind(find(ind>i12,1,'first'));
        str2 = text(i12:i2);
        str = sprintf('backgrounder%s',str2(11:end));
        eval(str)
        
        clear cStim 
%         backgrounder = background;
        for jj = 1:length(chUnitNum)
            cStim.name.Text = 'Spot Stimulus';
            cStim.cycles.Text = nCycles;
            cStim.control.wait{1}.Text = onTime+offTime;
            if length(chNum) == 1
                cStim.parameters.channel.Text= chNum;
            else
                cStim.parameters.channel.Text= chNum(jj);
            end
            cStim.parameters.channelUnit.Text= chUnitNum(jj);
            cStim.parameters.radius.Text= radius(jj);
            if length(brightness) == 1
                cStim.parameters.brightness.Text= brightness;
            else
                cStim.parameters.brightness.Text= brightness(jj);
            end
            if iscell(backgrounder)&length(backgrounder)>1
                cStim.parameters.background.Text= backgrounder{jj};
            else
                cStim.parameters.background.Text= backgrounder;
            end
            cStim.parameters.stimulusType.Text= char(stimType(jj));
            cStim.parameters.onTime.Text= onTime;
            cStim.parameters.offTime.Text= offTime;
            stimFile = [stimFile cStim];
            clear cStim
            cStim.name.Text = 'Spot Spontaneous';
            cStim.cycles.Text = nCycles;
            cStim.control.wait{1}.Text = spontaneous/nCycles;
            cStim.parameters.background.Text= backgrounder;
            stimFile = [stimFile cStim];
            if jj == length(chUnitNum)
                stimFile = [stimFile buffer];
            end
        end
        
        
    else
        cLogName = strcat(cFileName(1:end-4),'Log.xml');
        XMLfile = strcat(fN,cFileName);
        LogFile = strcat(fN,cLogName);
        if (exist(LogFile,'file')~=0&flag==1)|(exist(LogFile,'file')~=0&flag==2)
            logFileStruct = xml2struct(LogFile);
            lenLog = length(logFileStruct.logFile.Log);
            if lenLog > 1
                counter = 1;
                for jj = 1:lenLog
                    if strcmp(logFileStruct.logFile.Log{jj}.status.Text,'Normal')
                        endPoints(counter) = jj;
                        counter = counter+1;
                    end
                end
            else
                endPoints = 1;
            end

        else
            endPoints = 1;
        end
        cFile = xml2struct(XMLfile);
    
        for jj = endPoints
            if isfield(cFile.stimulusFile,'stimulus')
                if isempty(stimFile)
                    stimFile = cFile.stimulusFile.stimulus;
                else
                    stimFile = [stimFile cFile.stimulusFile.stimulus];
                end

                if ii ~= length(stimulusFiles)|jj~=endPoints(end)
                    stimFile = [stimFile buffer];
                end

            end
        end
    end   
end

if flag==1&sum(strcmp(stimulusFiles,'stimulusFile.xml'))==1
    cFileName = char(stimulusFiles(strcmp(stimulusFiles,'stimulusFile.xml')));
    XMLfile = strcat(fN,cFileName);
    filer = xml2struct(XMLfile);
    newFileName = strcat(fN,'stimulusFileOld.xml');
    struct2xml(filer,newFileName)
end



newFileName = strcat(fN,'stimulusFile.xml');
filer.stimulusFile.stimulus = stimFile;
struct2xml(filer,newFileName)



