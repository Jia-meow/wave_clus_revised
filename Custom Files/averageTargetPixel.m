function [ targetPixel,um,electrodes ] = averageTargetPixel(fileName)
% For LCD visual stimulation, this function calculates where the
% electrodes are located with respect to the LCD pixels

    electrodeList = [12:17, 21:28,31:38,41:48,51:58,61:68,71:78,82:87];
    calFile1 = loader(sprintf('%s/preCalDat.mat',fileName));
    calFile2 = loader(sprintf('%s/postCalDat.mat',fileName));
    if exist(calFile1) >0
        [targetPixel1,um1] = getPixelLocation(calFile1,electrodeList);
    else
        targetPixel1 = NaN(length(electrodeList),2);
        um1 = NaN;
    end
    if exist(calFile2) >0
        [targetPixel2,um2] = getPixelLocation(calFile2,electrodeList);
    else
        targetPixel2 = NaN(length(electrodeList),2);
        um2 = NaN;
    end
    if nanmax(nanmax(abs(targetPixel1-targetPixel2)))>=50
        targetPixel = targetPixel2;
        um = um2;
    else
        targetPixel(:,1) = nanmean([targetPixel1(:,1) targetPixel2(:,1)],2);
        targetPixel(:,2) = nanmean([targetPixel1(:,2) targetPixel2(:,2)],2);
        um = nanmean([um1 um2]);
    end
    if isnan(targetPixel(1))
%         error('No Calibration Data!')
    end
    targetPixel = targetPixel./25;
    electrodes = electrodeList;
end

