function preProcessData(fileName)
% Preprocesses the data by extracting the raw spike information and
% obtaining the rising and falling edges of the light and injection outputs
% as well as the stimulus tags (D1), cycle tags (D2), and pulse tags (D3)
% and saves these into the original file (fileName)


rawData = load(fileName);
% fields = fieldnames(rawData);
fieldNs = fields(rawData);
rawData = rmfield(rawData,fieldNs(strncmp(fieldNs,'adc',3)));
MCDfile = strcat(fileName(1:find(fileName == '.',1,'last')),'mcd');
[A1,A2,D1,D2,D3,time] = readMCD(MCDfile);
rawData.A1 = A1;
rawData.A2 = A2;
rawData.D1 = D1;
rawData.D2 = D2;
rawData.D3 = D3;
% rawData.time = time;
rawData.maxTime = max(time);
save(fileName,'-struct','rawData')


end

