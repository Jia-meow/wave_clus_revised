function [stimulusTypes,names] = readStimulusTypes(fileName)
% Read all available stimulus types


if nargin == 0
    fileName = loader('Labview/stimulusTypes.xml',1);
end

% fileName = 'stimulusTypes.xml';
xDoc = xmlread(fileName);
xRoot = xDoc.getDocumentElement;
stimuli = xRoot.getChildNodes;
stimNums = 1:2:(stimuli.getLength-2);
allNames = '';
for ii = 1:length(stimNums)
    firstStim = stimuli.item(stimNums(ii)).getChildNodes;
    node = firstStim.getFirstChild;
    numNodes = firstStim.getLength;
    counter = 1;
    name = {};
    data = {};
    for jj = 0:numNodes-1
        node = firstStim.item(jj);
        if strcmpi(node.getNodeName,'#text')

        else
            name(counter) = cellstr(char(node.getNodeName));
            data(counter) = cellstr(char(node.getTextContent));
            counter = counter + 1;
        end
    end
    nameInd = strcmp(name,'name');
    stimulusTypes(ii).name = char(data(nameInd));
    allNames = [allNames; cellstr(char(data(nameInd)))];
    commandInd = find(nameInd==0);
    for jj = 1:length(commandInd)
        stimulusTypes(ii).command(jj) = data(commandInd(jj));
    end
end

names = allNames;
end
