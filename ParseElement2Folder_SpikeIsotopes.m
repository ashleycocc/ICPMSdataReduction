%% Parse a folder of Element2 data files, perform on peak zeros subtraction
%  and calculate isotope ratios.
%  For 229Th-232Th-235U-238U runs

%% Setup for parsing text file
delimiter = '\t';
startRow = 6;
formatSpec = '%f%f%f%f%f%[^\n\r]';


%% Sort files out of folder

datafolderstring = './data/080118-AC-SpikeIsotopes/';
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

% remove 1st standard


%% Subtract on-peak zeroes
% assume all odd numbered rows are OPZs, even numbered rows are samples

isSTD = ~logical(1:length(data))';
is229 = ~logical(1:length(data))';
is235 = ~logical(1:length(data))';

for iRun = 2:2:length(data)
    
    % OPZ subtraction
    data(iRun).OPZcorrInt = data(iRun).intensities(50:end,:) - data(iRun-1).intensities(50:end,:);
    
    % Flag standards
    temp.nameString = data(iRun).fileName(4:6); 
    if strcmp(temp.nameString, 'std')
        isSTD(iRun) = 1;
    end
        
    %Flag Th-229 samples
    temp.sampNameString = data(iRun).fileName(4:6);
    if strcmp(temp.sampNameString, '229')
        is229(iRun) = 1;    
    end
    
    if strcmp(temp.sampNameString, '235')
        is235(iRun) = 1;    
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
SRMU970URatio = 0.00535;
mass238 = 238;
mass235 = 235;
stdData = data(isSTD);
Th229Data = data(is229);
U235Data = data(is235);

for iStd = 1:length(stdData)
    
    stdData(iStd).Beta = (log(CRM112aURatio) - log(stdData(iStd).U238U235Ratio)) ./ log(mass235/mass238);
    stdData(iStd).BetaAvg = mean(stdData(iStd).Beta);
    stdData(iStd).BetaStdError = std(stdData(iStd).Beta) / sqrt(length(stdData(iStd).Beta));
    stdData(iStd).MassFrCorrAvg238U235U = stdData(iStd).avg238U235URatios * (mass235/mass238) .^ stdData(iStd).BetaAvg;
    
    hold on
    plotHandle = errorbar(iStd, stdData(iStd).BetaAvg, stdData(iStd).BetaStdError, '.', 'MarkerEdgeColor','red','MarkerFaceColor','red');
    plotHandle.MarkerSize = 18;
    plotHandle.LineWidth = 2;
    plotHandle.Color = 'm';
    
end

for i235 = 1:length(U235Data)
    
    U235Data(i235).Beta = (log(SRMU970URatio) - log(U235Data(i235).U238U235Ratio)) ./ log(mass235/mass238);
    U235Data(i235).BetaAvg = mean(U235Data(i235).Beta);
    U235Data(i235).BetaStdError = std(U235Data(i235).Beta) / sqrt(length(U235Data(i235).Beta));
    U235Data(i235).MassFrCorrAvg238U235U = U235Data(i235).avg238U235URatios * (mass235/mass238) .^ U235Data(i235).BetaAvg;
    U235Data(i235).StdError238U235U = std(U235Data(i235).U238U235Ratio) / sqrt(length(U235Data(i235).U238U235Ratio));
    
    plotHandle2 = errorbar(i235, U235Data(i235).BetaAvg, U235Data(i235).BetaStdError, '.', 'MarkerEdgeColor','b','MarkerFaceColor','red');
    plotHandle2.MarkerSize = 18;
    plotHandle2.LineWidth = 2;
    plotHandle2.Color = 'c';
    
end


%% Apply Betas to spike cal data

samp12BetaAvg = mean([stdData(1:2).BetaAvg]);
samp34BetaAvg = mean([stdData(3:4).BetaAvg]);
samp56BetaAvg = mean([stdData(5:6).BetaAvg]);
samp7BetaAvg = mean([stdData(8:9).BetaAvg]);
mass229 = 229;
mass232 = 232;

for iMass12 = 1:2

    Th229Data(iMass12).massFracCorr232Th229Th = Th229Data(iMass12).Th232Th229Ratio * (mass229 / mass232 ) ^ samp12BetaAvg;
    Th229Data(iMass12).avgMassFracCorr232Th229Th = mean(Th229Data(iMass12).massFracCorr232Th229Th);
    Th229Data(iMass12).stdError232Th229Th = ((mean(std(Th229Data(iMass12).massFracCorr232Th229Th)...
        / sqrt(length(Th229Data(iMass12).massFracCorr232Th229Th)))));
    
    Th229Data(iMass12).massFracCorr238U235U = Th229Data(iMass12).U238U235Ratio * (mass235 / mass238 ) ^ samp12BetaAvg;
    Th229Data(iMass12).avgMassFracCorr238U235U = mean(Th229Data(iMass12).massFracCorr238U235U);
    
    U235Data(iMass12).massFracCorr232Th229Th = U235Data(iMass12).Th232Th229Ratio * (mass229 / mass232 ) ^ samp12BetaAvg;
    U235Data(iMass12).avgMassFracCorr232Th229Th = mean(U235Data(iMass12).massFracCorr232Th229Th);
    
    U235Data(iMass12).massFracCorr238U235U = U235Data(iMass12).U238U235Ratio * (mass235 / mass238 ) ^ samp12BetaAvg;
    U235Data(iMass12).avgMassFracCorr238U235U = mean(U235Data(iMass12).massFracCorr238U235U);
    
end


for iMass34 = 3:4

    Th229Data(iMass34).massFracCorr232Th229Th = Th229Data(iMass34).Th232Th229Ratio * (mass229 / mass232 ) ^ samp34BetaAvg;
    Th229Data(iMass34).avgMassFracCorr232Th229Th = mean(Th229Data(iMass34).massFracCorr232Th229Th);
    Th229Data(iMass34).stdError232Th229Th = (((mean(std(Th229Data(iMass34).massFracCorr232Th229Th)...
        / sqrt(length(Th229Data(iMass34).massFracCorr232Th229Th)) ))));
    
    Th229Data(iMass34).massFracCorr238U235U = Th229Data(iMass34).U238U235Ratio * (mass235 / mass238 ) ^ samp34BetaAvg;
    Th229Data(iMass34).avgMassFracCorr238U235U = mean(Th229Data(iMass34).massFracCorr238U235U);
    
    U235Data(iMass34).massFracCorr232Th229Th = U235Data(iMass34).Th232Th229Ratio * (mass229 / mass232 ) ^ samp12BetaAvg;
    U235Data(iMass34).avgMassFracCorr232Th229Th = mean(U235Data(iMass34).massFracCorr232Th229Th);
    
    U235Data(iMass34).massFracCorr238U235U = U235Data(iMass34).U238U235Ratio * (mass235 / mass238 ) ^ samp12BetaAvg;
    U235Data(iMass34).avgMassFracCorr238U235U = mean(U235Data(iMass34).massFracCorr238U235U);
    
end

for iMass56 = 5:6

    Th229Data(iMass56).massFracCorr232Th229Th = Th229Data(iMass56).Th232Th229Ratio * (mass229 / mass232 ) ^ samp56BetaAvg;
    Th229Data(iMass56).avgMassFracCorr232Th229Th = mean(Th229Data(iMass56).massFracCorr232Th229Th);
    Th229Data(iMass56).stdError232Th229Th = (((mean(std(Th229Data(iMass56).massFracCorr232Th229Th)...
        / sqrt(length(Th229Data(iMass56).massFracCorr232Th229Th)) ))));
    
    Th229Data(iMass56).massFracCorr238U235U = Th229Data(iMass56).U238U235Ratio * (mass235 / mass238 ) ^ samp56BetaAvg;
    Th229Data(iMass56).avgMassFracCorr238U235U = mean(Th229Data(iMass56).massFracCorr238U235U);
    
    U235Data(iMass56).massFracCorr232Th229Th = U235Data(iMass56).Th232Th229Ratio * (mass229 / mass232 ) ^ samp12BetaAvg;
    U235Data(iMass56).avgMassFracCorr232Th229Th = mean(U235Data(iMass56).massFracCorr232Th229Th);
    
    U235Data(iMass56).massFracCorr238U235U = U235Data(iMass56).U238U235Ratio * (mass235 / mass238 ) ^ samp12BetaAvg;
    U235Data(iMass56).avgMassFracCorr238U235U = mean(U235Data(iMass56).massFracCorr238U235U);
    
end

for iMass7 = 7

    Th229Data(iMass7).massFracCorr232Th229Th = Th229Data(iMass7).Th232Th229Ratio * (mass229 / mass232 ) ^ samp7BetaAvg;
    Th229Data(iMass7).avgMassFracCorr232Th229Th = mean(Th229Data(iMass7).massFracCorr232Th229Th);
    Th229Data(iMass7).stdError232Th229Th = (((mean(std(Th229Data(iMass7).massFracCorr232Th229Th)...
        / sqrt(length(Th229Data(iMass7).massFracCorr232Th229Th)) ))));
    
    Th229Data(iMass7).massFracCorr238U235U = Th229Data(iMass7).U238U235Ratio * (mass235 / mass238 ) ^ samp7BetaAvg;
    Th229Data(iMass7).avgMassFracCorr238U235U = mean(Th229Data(iMass7).massFracCorr238U235U);
    
    U235Data(iMass7).massFracCorr232Th229Th = U235Data(iMass7).Th232Th229Ratio * (mass229 / mass232 ) ^ samp12BetaAvg;
    U235Data(iMass7).avgMassFracCorr232Th229Th = mean(U235Data(iMass7).massFracCorr232Th229Th);
    
    U235Data(iMass7).massFracCorr238U235U = U235Data(iMass7).U238U235Ratio * (mass235 / mass238 ) ^ samp12BetaAvg;
    U235Data(iMass7).avgMassFracCorr238U235U = mean(U235Data(iMass7).massFracCorr238U235U);
    
end


mean232229Thratio = mean(Th229Data.massFracCorr232Th229Th);


