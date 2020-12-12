%% Parse a folder of Element2 data files, perform on peak zeros subtraction
%  and calculate isotope ratios.
%  For UTEVA column calibration

%% Setup for parsing text file
delimiter = '\t';
startRow = 6;
formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';


%% Sort files out of folder

datafolderstring = './data/080518-AC-UTEVA-Calib-2/';
% don't forget to put a / at the end of the folder name!

fileStruct = dir([datafolderstring '*.txt']);
n.files = numel(fileStruct) - 1;
% minus one for the .txt file with the sample list generated from sequence
[temp.min, temp.indx] = min([fileStruct.bytes]); % sample list is smallest file
fileStruct(temp.indx) = []; % delete sample list

data = struct('fileName', {}, 'timeStamp', {}, 'intensities', {});


%% Pull data out of files

for iFile = 1:n.files
    
    % extract filename
    temp.fileName = fileStruct(iFile).name;
    data(iFile).fileName = temp.fileName(1:end-4);
    
    % extract timestamp as MATLAB time
    data(iFile).timeStamp = fileStruct(iFile).datenum;
    
    % extract intensities as four-column matrix
    % where columns are 229Th - 232Th - 235U - 238U
    fileID = fopen([datafolderstring fileStruct(iFile).name]);
    
   textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', ...
        'ReturnOnError', false, 'EndOfLine', '\r\n');
    temp.dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
        'TextType', 'string', 'ReturnOnError', false);
    
    data(iFile).intensities =  [temp.dataArray{1:end-1}];
    data(iFile).intensities = data(iFile).intensities(:,2:end);
    
    fclose(fileID);
    
end

% sort by timestamp
[~,temp.sortIndex] = sortrows([data.timeStamp].'); 
data = data(temp.sortIndex);

%% Subtract on-peak zeroes
% assume all odd numbered rows are OPZs, even numbered rows are samples

isSample = ~logical(1:length(data))';
isUstd = ~logical(1:length(data))';

for iRun = 2:2:length(data)
    
    % OPZ subtraction
    data(iRun).OPZcorrInt = data(iRun).intensities(20:81,:) - data(iRun-1).intensities(20:81,:);

        
    %Flag U stds 
    temp.nameString = data(iRun).fileName(4:5);
    if strcmp(temp.nameString, '99')
        isUstd(iRun) = 1;
    end
    
    %Flag Samples
    temp.sampNameString = data(iRun).fileName(10:end);
    if strcmp(temp.sampNameString, 'samp')
        isSample(iRun) = 1;
    elseif strcmp(temp.sampNameString, 'conv')
        isSample(iRun) = 1;
    elseif strcmp(temp.sampNameString, 'U')
        isSample(iRun) = 1;    
    elseif strcmp(temp.sampNameString, 'Th')
        isSample(iRun) = 1;
    elseif strcmp(temp.sampNameString, 'rinse')
        isSample(iRun) = 1;
   
    end    
    
    
end


%% calculate average intensities


sampData = data(isSample);
stdData = data(isUstd);
  
for iSample= 1:length(sampData)
    
    sampData(iSample).fileName = sampData(iSample).fileName(4:end);
    
end

 [~,temp.sampNameIndex] = sortrows({sampData.fileName}.'); 
 sampData = sampData(temp.sampNameIndex); clear index

ColVol = "UTEVAColVol.xlsx";
sheet = 1;
CVexcel = 'D5:D25';
CV = xlsread(ColVol, CVexcel);

 for iCV = 1:length(sampData)
 
     sampData(iCV).ColumnVolume = CV;
     
 end   
 
 for iCVind = 1:length(sampData)
 
     sampData(iCVind).ColumnVolumeInd = sampData(iCVind).ColumnVolume(iCVind, 1);
     
 end
 
 for iCVWash = 1:11
 
     sampData(iCVWash).ColumnVolumeWash = sampData(iCVWash).ColumnVolume(iCVWash, 1);
 
 end
     
 for iEluteSampRatio = 12:length(sampData)
     
     sampData(iEluteSampRatio).Ca = (sampData(iEluteSampRatio).OPZcorrInt(:,2));
     sampData(iEluteSampRatio).Ti = (sampData(iEluteSampRatio).OPZcorrInt(:,3));
     sampData(iEluteSampRatio).Th232 = (sampData(iEluteSampRatio).OPZcorrInt(:,5));
     sampData(iEluteSampRatio).U238 = (sampData(iEluteSampRatio).OPZcorrInt(:,6));
     
     sampData(iEluteSampRatio).avgCa = mean(sampData(iEluteSampRatio).Ca);
     sampData(iEluteSampRatio).avgTi = mean(sampData(iEluteSampRatio).Ti);
     sampData(iEluteSampRatio).avg238U = mean(sampData(iEluteSampRatio).U238);
     sampData(iEluteSampRatio).avg232Th = mean(sampData(iEluteSampRatio).Th232) ;
     
 end
 
 for iWashSampRatio = 1:11
     
     sampData(iWashSampRatio).Mg = (sampData(iWashSampRatio).OPZcorrInt(:,1));
     sampData(iWashSampRatio).Ca = (sampData(iWashSampRatio).OPZcorrInt(:,2));
     sampData(iWashSampRatio).Ti = (sampData(iWashSampRatio).OPZcorrInt(:,3));
     sampData(iWashSampRatio).Zr = (sampData(iWashSampRatio).OPZcorrInt(:,5));
     sampData(iWashSampRatio).Th232 = (sampData(iWashSampRatio).OPZcorrInt(:,6));
     sampData(iWashSampRatio).U238 = (sampData(iWashSampRatio).OPZcorrInt(:,7));
     
     
     sampData(iWashSampRatio).avgMg = mean(sampData(iWashSampRatio).Mg);
     sampData(iWashSampRatio).avgCa = mean(sampData(iWashSampRatio).Ca);
     sampData(iWashSampRatio).avgTi = mean(sampData(iWashSampRatio).Ti);
     sampData(iWashSampRatio).avgZr = mean(sampData(iWashSampRatio).Zr);
     sampData(iWashSampRatio).avg238U = mean(sampData(iWashSampRatio).U238) ;
     sampData(iWashSampRatio).avg232Th = mean(sampData(iWashSampRatio).Th232) ;
     
 end
 
 sumMgintensity = sum([sampData(:).avgMg]);
 sumCaintensity = sum([sampData(:).avgCa]);
 sumTiintensity = sum([sampData(:).avgTi]);
 sumZrintensity = sum([sampData(:).avgZr]);
 sum232Thintensity = sum([sampData(:).avg232Th]);
 sum238Uintensity = sum([sampData(:).avg238U]);
 
 for iYield = 1:length(sampData)
    
    sampData(iYield).PercentMgYield = (sampData(iYield).avgMg / sumMgintensity) * 100;
    sampData(iYield).PercentCaYield = (sampData(iYield).avgCa / sumCaintensity) * 100;
    sampData(iYield).PercentTiYield = (sampData(iYield).avgTi / sumTiintensity) * 100;
    sampData(iYield).PercentZrYield = (sampData(iYield).avgZr / sumZrintensity) * 100; 
    sampData(iYield).Percent232ThYield = (sampData(iYield).avg232Th / sum232Thintensity) * 100;
    sampData(iYield).Percent238UYield = (sampData(iYield).avg238U / sum238Uintensity) * 100;
    
 end


 
%% Create Column Calibration curve

hold on
plotHandle = plot([sampData.ColumnVolumeInd], [sampData.Percent238UYield]);
plotHandle.LineWidth = 2;
% plotHandle.Color = 'b';

plotHandle2 = plot([sampData.ColumnVolumeInd], [sampData.Percent232ThYield]);
plotHandle2.LineWidth = 2;
% plotHandle2.Color = 'c';

plotHandle3 = plot([sampData.ColumnVolumeWash], [sampData.PercentZrYield]);
plotHandle3.LineWidth = 2;
% plotHandle3.Color = 'm';

plotHandle4 = plot([sampData.ColumnVolumeInd], [sampData.PercentTiYield]);
plotHandle4.LineWidth = 2;
% plotHandle4.Color = 'r';

plotHandle5 = plot([sampData.ColumnVolumeInd], [sampData.PercentCaYield]);
plotHandle5.LineWidth = 2;
% plotHandle5.Color = 'y';

plotHandle6 = plot([sampData.ColumnVolumeWash], [sampData.PercentMgYield]);
plotHandle6.LineWidth = 2;
% plotHandle6.Color = 'k';


ax = gca;
ax.XLabel.String = 'Aliquot (uL)';
ax.XLabel.FontSize = 18;
ax.YLabel.String = 'Percent Yield';
ax.YLabel.FontSize = 18;
ax.XLim = [0 3050];
ax.XTick = [0:125:3050];


