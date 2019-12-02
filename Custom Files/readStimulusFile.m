function stimulusFile = readStimulusFile(fileName)

if ~strcmp(fileName(end-2:end),'xml')
    fileName = loader(sprintf('%s/stimulusFile.xml',fileName));
end

% fileName = 'stimulusFile.xml';
xDoc = xmlread(fileName);
xRoot = xDoc.getDocumentElement;
stimuli = xRoot.getChildNodes;
stimNums = 1:2:(stimuli.getLength-2);
for ii = 1:length(stimNums)
    currentStim = stimuli.item(stimNums(ii)).getChildNodes;
    numNodes = currentStim.getLength;
    name = '';
    cycleNum = 1;
    control = '';
    controlVal = '';
    parameter = '';
    parameterVal = '';
    for jj = 0:numNodes-1
        node = currentStim.item(jj);
        if strcmpi(node.getNodeName,'#text')
            
        else
            switch (char(node.getNodeName))
                case 'name'
                    name = char(node.getTextContent);
                case 'cycles'
                    cycleNum = str2num(node.getTextContent);
                case 'control'
                    controlNodes = node.getChildNodes;
                    for kk = 0:controlNodes.getLength-1
                        if ~strncmp(controlNodes.item(kk).getNodeName,'#',1)
                            control = [control; cellstr(char(controlNodes.item(kk).getNodeName))];
                            controlVal = [controlVal; cellstr(char(controlNodes.item(kk).getTextContent))];
                        end
                    end
                case 'parameters'
                    parameterNodes = node.getChildNodes;
                    for kk = 0:parameterNodes.getLength-1
                        if ~strncmp(parameterNodes.item(kk).getNodeName,'#',1)
                            parameter = [parameter; cellstr(char(parameterNodes.item(kk).getNodeName))];
                            parameterVal = [parameterVal; cellstr(char(parameterNodes.item(kk).getTextContent))];
                        end
                    end                    
            end
        end
    end
    stimulusFile(ii).name = name;
    stimulusFile(ii).cycleNum = cycleNum;
    stimulusFile(ii).controlNames = control;
    stimulusFile(ii).controlVals = controlVal;
    stimulusFile(ii).parameterNames = parameter;
    stimulusFile(ii).parameterVals = parameterVal;
end
if isempty(stimNums)
    stimulusFile = NaN;
end
end