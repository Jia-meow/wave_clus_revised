function [varargout] = getUnitFiles(folderName,dirFlag)
%Searches directory and outputs all of the units, unit file names,
%channels, and channel units

if nargin == 1
    dirFlag = true;
end

if iscell(folderName)
    folderName = char(folderName);
end


if length(strfind(folderName,'/'))==1 & dirFlag == true
    folderName = loader(sprintf('%s/processedData',folderName));
elseif dirFlag == false
    folderName = loader(sprintf('%s',folderName)); 
else
    folderName = sprintf('%s/processedData',folderName);
end


unitDir = dir(folderName);
for ii = 1:length(unitDir)
    unitFiles{ii} = unitDir(ii).name;
end
unitFiles = unitFiles(strncmp(unitFiles,'u',1));

% if dirFlag == true
    unitFiles = arrayfun(@(col) sprintf('%s/%s', folderName, char(unitFiles{col})), 1:size(unitFiles,2),'Unif', false)';
    
varargout{1} = unitFiles;
if nargout == 2
    strSplitter= regexp(unitFiles,'unit_','split');
    units = cell2mat(arrayfun(@(x) str2num(strSplitter{x}{2}(1:3)),1:length(unitFiles),'Unif',false));
    varargout{2} = units;
end

if nargout == 3
    strSplitter= regexp(unitFiles,'unit_','split');
    units = cell2mat(arrayfun(@(x) str2num(strSplitter{x}{2}(1:3)),1:length(unitFiles),'Unif',false));
    varargout{2} = units;
    for ii = 1:length(unitFiles)
        dat = load(char(unitFiles(ii)));
        chNum(ii) = dat.unit.channelNumber;
        
    end
    varargout{3} = chNum;
end

if nargout == 4
    strSplitter= regexp(unitFiles,'unit_','split');
    units = cell2mat(arrayfun(@(x) str2num(strSplitter{x}{2}(1:3)),1:length(unitFiles),'Unif',false));
    varargout{2} = units;
    for ii = 1:length(unitFiles)
        dat = load(char(unitFiles(ii)));
        chNum(ii) = dat.unit.channelNumber;
        fs = fields(dat.unit);
        if sum(strcmp(fs,'channelUnitNumber'))~=0
            chUnit(ii) = dat.unit.channelUnitNumber;
        else
            chUnit(ii) = NaN;
        end
    end
    varargout{3} = chNum;
    varargout{4} = chUnit;
end

% end

