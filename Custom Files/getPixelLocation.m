function [varargout] = getPixelLocation(calFile,ChOI)
%Find which pixel of an LCD screen that the current channel corresponds to

if nargin == 1
    ChOI = [12:17, 21:28,31:38,41:48,51:58,61:68,71:78,82:87];
end

if ~strcmp(calFile(end-3:end),'.mat')
    folderName = calFile;
    calFile1 = loader(sprintf('%s/preCalDat.mat',folderName));
    calFile2 = loader(sprintf('%s/postCalDat.mat',folderName));
    if exist(calFile1) >0
        [targetPixel1,um1] = getPixelLocation(calFile1,ChOI);
    else
        targetPixel1 = NaN(length(ChOI),2);
        um1 = NaN;
    end
    if exist(calFile2) >0
        [targetPixel2,um2] = getPixelLocation(calFile2,ChOI);
    else
        targetPixel2 = NaN(length(ChOI),2);
        um2 = NaN;
    end
    if nanmax(nanmax(abs(targetPixel1-targetPixel2)))>=50
        targetPixel = targetPixel2;
        umPixel = um2;
    else
        targetPixel(:,1) = nanmean([targetPixel1(:,1) targetPixel2(:,1)],2);
        targetPixel(:,2) = nanmean([targetPixel1(:,2) targetPixel2(:,2)],2);
        umPixel = nanmean([um1 um2]);
    end
    if isnan(targetPixel(1))
        error('No Calibration Data!')
    end
    targetPixel = targetPixel./25;

    
else


    electrodeSize = 30;
    interElectrodeDistance = 200;

    % calDat = load(loader('airPulseTimingTest/exCalDataur.mat'));
    % calDat = load(loader('20140528/Recording1/preCalDat.mat'));
    calDat = load(calFile);

    calDat = calDat.calDat;


    fieldNames = fields(calDat);
    if sum(strcmp(fieldNames,'note'))>0
       elecInds = strncmp(fieldNames,'elec',4);
       electrodeNames = fieldNames(elecInds);
       for ii = 1:length(electrodeNames)
           electrodeNumbers(ii) = str2num(electrodeNames{ii}(end-1:end));
       end
       notes = calDat.note;
       for ii = 1:length(notes(:,1))
           ind = find(electrodeNumbers==notes(ii,1));
           electrodeNumbers(ind) = notes(ii,2);
       end
    else
       elecInds = strncmp(fieldNames,'elec',4);
       electrodeNames = fieldNames(elecInds);
       for ii = 1:length(electrodeNames)
           electrodeNumbers(ii) = str2num(electrodeNames{ii}(end-1:end));
       end
    end

    [X,Y] = meshgrid(1:size(calDat.electrode21,2),1:size(calDat.electrode21,1));

    electrodes = 11:88;
    electrodeLocs = nan(length(electrodes),2);
    electrodeLocs(electrodes==21,:) = [mean(X((calDat.electrode21>1))),mean(Y((calDat.electrode21>1)))];
    electrodeLocs(electrodes==28,:) = [mean(X((calDat.electrode28>1))),mean(Y((calDat.electrode28>1)))];
    electrodeLocs(electrodes==71,:) = [mean(X((calDat.electrode71>1))),mean(Y((calDat.electrode71>1)))];
    electrodeLocs(electrodes==78,:) = [mean(X((calDat.electrode78>1))),mean(Y((calDat.electrode78>1)))];

    clear len1 len2 
    len1(1,:) = [electrodeLocs(electrodes==28,1)-electrodeLocs(electrodes==21,1) , electrodeLocs(electrodes==28,2)-electrodeLocs(electrodes==21,2)];
    len1(2,:) = [electrodeLocs(electrodes==78,1)-electrodeLocs(electrodes==71,1) , electrodeLocs(electrodes==78,2)-electrodeLocs(electrodes==71,2)];
    len2(1,:) = [electrodeLocs(electrodes==71,1)-electrodeLocs(electrodes==21,1) , electrodeLocs(electrodes==71,2)-electrodeLocs(electrodes==21,2)];
    len2(2,:) = [electrodeLocs(electrodes==78,1)-electrodeLocs(electrodes==28,1) , electrodeLocs(electrodes==78,2)-electrodeLocs(electrodes==28,2)];

    e1X = floor(electrodeNumbers(1)/10);
    e1Y = electrodeNumbers(1)-e1X*10;
    e2X = floor(electrodeNumbers(2)/10);
    e2Y = electrodeNumbers(2)-e2X*10;
    e3X = floor(electrodeNumbers(3)/10);
    e3Y = electrodeNumbers(3)-e3X*10;
    e4X = floor(electrodeNumbers(4)/10);
    e4Y = electrodeNumbers(4)-e4X*10;

    refDist1(1) = sqrt( ((e1Y-e3Y)*interElectrodeDistance).^2 + ((e1X-e3X)*interElectrodeDistance).^2);
    refDist1(2) = sqrt( ((e2Y-e4Y)*interElectrodeDistance).^2 + ((e2X-e4X)*interElectrodeDistance).^2);
    refDist2(1) = sqrt( ((e1Y-e2Y)*interElectrodeDistance).^2 + ((e1X-e2X)*interElectrodeDistance).^2);
    refDist2(2) = sqrt( ((e3Y-e4Y)*interElectrodeDistance).^2 + ((e3X-e4X)*interElectrodeDistance).^2);


    dist1 = sqrt(len1(:,1).^2+len1(:,2).^2)';
    dist2 = sqrt(len2(:,1).^2+len2(:,2).^2)';

    umPixel = mean(mean([refDist1./dist1; refDist2./dist2]));

    refComp(1,:) =  [((e3X-e1X)*interElectrodeDistance) ((e3Y-e1Y)*interElectrodeDistance)];
    refComp(2,:) =  [((e4X-e2X)*interElectrodeDistance) ((e4Y-e2Y)*interElectrodeDistance)];
    refComp(3,:) =  [((e2X-e1X)*interElectrodeDistance) ((e2Y-e1Y)*interElectrodeDistance)];
    refComp(4,:) =  [((e4X-e3X)*interElectrodeDistance) ((e4Y-e3Y)*interElectrodeDistance)];
    pixComp = [len1;len2];

    pixAngle = rad2deg(atan2(pixComp(:,2),pixComp(:,1)));
    refAngle = rad2deg(atan2(refComp(:,2),refComp(:,1)));
    diffAngle = mean(refAngle-pixAngle);


    for ii = 1:length(ChOI)
%         refCh = 21;
        % targCh = 28;
        refCh = electrodeNumbers(1);
        targCh = ChOI(ii);

        rX = floor(refCh/10);
        rY = refCh-rX*10;
        tX = floor(targCh/10);
        tY = targCh-tX*10; 

        dX = tX-rX;
        dY = tY-rY;
        refDist = sqrt((dX*interElectrodeDistance).^2 + (dY*interElectrodeDistance).^2);
        refLoc = electrodeLocs(electrodes==21,:)*umPixel;
        refAngle = rad2deg(atan2(dY,dX))-diffAngle;

        targetLoc = refLoc+[refDist*cos(deg2rad(refAngle)) refDist*sin(deg2rad(refAngle))];
        targetPixel(ii,:) = targetLoc./umPixel;
    end

end
varargout{1} = targetPixel;
if nargout == 2
    varargout{2} = umPixel;
elseif nargout == 3
    varargout{2} = umPixel;
    varargout{3} = ChOI;
end












end

