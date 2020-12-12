function sampData = ParseElement2Folder_SpikeCalFunc;

%% Parse a folder of Element2 data files, perform on peak zeros subtraction
%  and calculate isotope ratios.
%  For 229Th-232Th-235U-238U runs

%% Setup for parsing text file
delimiter = '\t';
startRow = 6;
formatSpec = '%f%f%f%f%f%[^\n\r]';


%% Sort files out of folder

datafolderstring = './data/08252018_AC-SpikeCal/';
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

isSTD = ~logical(1:length(data))';
isSamp = ~logical(1:length(data))';

for iRun = 2:2:length(data)
    
    % OPZ subtraction
    data(iRun).OPZcorrInt = data(iRun).intensities(30:(end-20),:) - data(iRun-1).intensities(30:(end-20),:);
    
    % Flag standards
    temp.nameString = data(iRun).fileName(4:5); 
    if strcmp(temp.nameString, '99')
        isSTD(iRun) = 1;
    end
        
    %Flag samples
    temp.sampNameString = data(iRun).fileName(10:11);
    if strcmp(temp.sampNameString, 'SC')
        isSamp(iRun) = 1;
    end
    
     
    end
    


%% calculate ratios and average ratios


for iRatio = 2:2:length(data)

    Th232 = (data(iRatio).OPZcorrInt(:,2));
    Th229 = (data(iRatio).OPZcorrInt(:,1));
    data(iRatio).Th232Th229Ratio = Th232 ./ Th229 ;
    
    
    U235 = (data(iRatio).OPZcorrInt(:,3));
    U238 = (data(iRatio).OPZcorrInt(:,4));
    data(iRatio).U238U235Ratio = U238 ./ U235 ;
    
    
    data(iRatio).avg232Th229ThRatios = mean(data(iRatio).Th232Th229Ratio);
    data(iRatio).avg238U235URatios = mean(data(iRatio).U238U235Ratio);
end

%% Calculate betas and correct average ratios for mass fractionation

CRM112aURatio = 137.844;
SRMU970URatio = 186.78;
mass238 = 238;
mass235 = 235;
stdData = data(isSTD);
sampData = data(isSamp);

for iStd = 1:length(stdData)
    
    stdData(iStd).Beta = (log(CRM112aURatio) - log(stdData(iStd).U238U235Ratio)) ./ log(mass235/mass238);
    stdData(iStd).BetaAvg = mean(stdData(iStd).Beta);
    stdData(iStd).BetaStdError = std(stdData(iStd).Beta) / sqrt(length(stdData(iStd).Beta));
    stdData(iStd).MassFrCorr238U235U = stdData(iStd).U238U235Ratio .* (mass235/mass238) .^ stdData(iStd).Beta;
    stdData(iStd).MassFrCorrAvg238U235U = mean(stdData(iStd).MassFrCorr238U235U);
    stdData(iStd).U238U235StdError = std(stdData(iStd).U238U235Ratio) / (sqrt(length(stdData(iStd).U238U235Ratio)));
    stdData(iStd).U238U235FracUncert = mean(stdData(iStd).U238U235StdError ./ stdData(iStd).U238U235Ratio);
    
%     hold on
%     plotHandle = errorbar(iStd, stdData(iStd).BetaAvg, stdData(iStd).BetaStdError, '.', 'MarkerEdgeColor','red','MarkerFaceColor','red');
%     plotHandle.MarkerSize = 18;
%     plotHandle.LineWidth = 2;
%     plotHandle.Color = 'm';
    
end

% hold off

%% Apply Betas to spike cal data

SC05xBetaAvg = mean([stdData(1:3).BetaAvg]);
SC1xBetaAvg = mean([stdData(5:7).BetaAvg]);
mass229 = 229;
mass232 = 232;

for iMass = 1:5

    sampData(iMass).massFracCorr232Th229Th = sampData(iMass).Th232Th229Ratio * (mass229 / mass232 ) ^ SC05xBetaAvg;
    sampData(iMass).avgMassFracCorr232Th229Th = mean(sampData(iMass).massFracCorr232Th229Th);
    sampData(iMass).stdError232Th229Th = (sqrt((mean(std(sampData(iMass).massFracCorr232Th229Th)...
        / sqrt(length(sampData(iMass).massFracCorr232Th229Th)) ./ sampData(iMass).massFracCorr232Th229Th))^2 ...
        + (stdData(iMass).U238U235FracUncert^2))) * sampData(iMass).avgMassFracCorr232Th229Th;
    
    sampData(iMass).massFracCorr238U235U = sampData(iMass).U238U235Ratio * (mass235 / mass238 ) ^ SC05xBetaAvg;
    sampData(iMass).avgMassFracCorr238U235U = mean(sampData(iMass).massFracCorr238U235U);
    sampData(iMass).stdError238U235U = (sqrt((mean(std(sampData(iMass).massFracCorr238U235U)...
        / sqrt(length(sampData(iMass).massFracCorr238U235U)) ./ sampData(iMass).massFracCorr238U235U))^2 ...
        + (stdData(iMass).U238U235FracUncert^2))) * sampData(iMass).avgMassFracCorr238U235U;
    
end


for iMass1x = 6:10

    sampData(iMass1x).massFracCorr232Th229Th = sampData(iMass1x).Th232Th229Ratio * (mass229 / mass232 ) ^ SC05xBetaAvg;
    sampData(iMass1x).avgMassFracCorr232Th229Th = mean(sampData(iMass1x).massFracCorr232Th229Th);
    sampData(iMass1x).stdError232Th229Th = (sqrt((mean(std(sampData(iMass1x).massFracCorr232Th229Th)...
        / sqrt(length(sampData(iMass1x).massFracCorr232Th229Th)) ./ sampData(iMass1x).massFracCorr232Th229Th))^2 ...
        + (stdData(iMass1x-2).U238U235FracUncert^2))) * sampData(iMass1x).avgMassFracCorr232Th229Th;
    
    sampData(iMass1x).massFracCorr238U235U = sampData(iMass1x).U238U235Ratio * (mass235 / mass238 ) ^ SC05xBetaAvg;
    sampData(iMass1x).avgMassFracCorr238U235U = mean(sampData(iMass1x).massFracCorr238U235U);
    sampData(iMass1x).stdError238U235U = (sqrt((mean(std(sampData(iMass1x).massFracCorr238U235U)...
        / sqrt(length(sampData(iMass1x).massFracCorr238U235U)) ./ sampData(iMass1x).massFracCorr238U235U))^2 ...
        + (stdData(iMass1x-2).U238U235FracUncert^2))) * sampData(iMass1x).avgMassFracCorr238U235U;
    
    
end



