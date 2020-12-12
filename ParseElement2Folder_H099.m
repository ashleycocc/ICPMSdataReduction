%% Parse a folder of Element2 data files, perform on peak zeros subtraction
%  and calculate isotope ratios.
%  For 229Th-232Th-235U-238U runs

%% Setup for parsing text file
delimiter = '\t';
startRow = 6;
formatSpec = '%f%f%f%f%f%[^\n\r]';


%% Sort files out of folder

datafolderstring = './data/022219_AC_H099_2/';
% don't forget to put a / at the end of the folder name!

fileStruct = dir([datafolderstring '*.txt']);
n.files = numel(fileStruct) -1;
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
isNotRinse = ~logical(1:length(data))';

for iRinse = 1:99
    
    % Remove rinses
    temp.nameString1 = data(iRinse).fileName(4:5); 
    if strcmp(temp.nameString1, '99')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '11')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '00')
        isNotRinse(iRinse) = 1;
    end

end

for iRinse2 = 100:length(data)
    
    % Remove rinses
    temp.nameString2 = data(iRinse2).fileName(5:6); 
    if strcmp(temp.nameString2, '99')
        isNotRinse(iRinse2) = 1;
    end
    if strcmp(temp.nameString2, '11')
        isNotRinse(iRinse2) = 1;
    end
    if strcmp(temp.nameString2, '00')
        isNotRinse(iRinse2) = 1;
    end

end

data = data(isNotRinse);

for iRun = 2:2:66
    
    % OPZ subtraction
    data(iRun).OPZcorrInt = data(iRun).intensities(30:(end-30),:) - data(iRun-1).intensities(30:(end-30),:);
    
    % Flag standards
    temp.nameString = data(iRun).fileName(4:5); 
    if strcmp(temp.nameString, '99')
        isSTD(iRun) = 1;
    end
        
    %Flag samples
    temp.sampNameString = data(iRun).fileName(4:5);
    if strcmp(temp.sampNameString, '11')
        isSamp(iRun) = 1;
    end
    
     
    end
    
for iRun2 = 66:2:length(data)
    
    % OPZ subtraction
    data(iRun2).OPZcorrInt = data(iRun2).intensities(30:(end-30),:) - data(iRun2-1).intensities(30:(end-30),:);
    
    % Flag standards
    temp.nameString3 = data(iRun2).fileName(5:6); 
    if strcmp(temp.nameString3, '99')
        isSTD(iRun2) = 1;
    end
        
    %Flag samples
    temp.sampNameString3 = data(iRun2).fileName(5:6);
    if strcmp(temp.sampNameString3, '11')
        isSamp(iRun2) = 1;
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
    
    data(iRatio).U235Mean = mean(U235);
    data(iRatio).Th229Mean = mean(Th229);
    data(iRatio).U238Mean = mean(U238);
    data(iRatio).Th232Mean = mean(Th232);
    
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

%% Apply Betas to data

SC05xBetaAvg = mean([stdData(1:6).BetaAvg]);
SC1xBetaAvg = mean([stdData(6:12).BetaAvg]);
mass229 = 229;
mass232 = 232;

for iMass = 1:18

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


for iMass = 19:34

    sampData(iMass).massFracCorr232Th229Th = sampData(iMass).Th232Th229Ratio * (mass229 / mass232 ) ^ SC05xBetaAvg;
    sampData(iMass).avgMassFracCorr232Th229Th = mean(sampData(iMass).massFracCorr232Th229Th);
    sampData(iMass).stdError232Th229Th = (sqrt((mean(std(sampData(iMass).massFracCorr232Th229Th)...
        / sqrt(length(sampData(iMass).massFracCorr232Th229Th)) ./ sampData(iMass).massFracCorr232Th229Th))^2 ...
        + (stdData(iMass-16).U238U235FracUncert^2))) * sampData(iMass).avgMassFracCorr232Th229Th;
    
    sampData(iMass).massFracCorr238U235U = sampData(iMass).U238U235Ratio * (mass235 / mass238 ) ^ SC05xBetaAvg;
    sampData(iMass).avgMassFracCorr238U235U = mean(sampData(iMass).massFracCorr238U235U);
    sampData(iMass).stdError238U235U = (sqrt((mean(std(sampData(iMass).massFracCorr238U235U)...
        / sqrt(length(sampData(iMass).massFracCorr238U235U)) ./ sampData(iMass).massFracCorr238U235U))^2 ...
        + (stdData(iMass-16).U238U235FracUncert^2))) * sampData(iMass).avgMassFracCorr238U235U;
    

end


%% Plotting ratios to see if they make sense

figHandle = figure('Name', '238U235URatios');
pH2 = subplot(5, 2, 1);
plotHandle2 = plot(sampData(1).massFracCorr238U235U, '.');
plotHandle2.MarkerSize = 14;
pH3 = subplot(5, 2, 2);
plotHandle3 = plot(sampData(3).massFracCorr238U235U, '.');
plotHandle3.MarkerSize = 14;
pH4 = subplot(5, 2, 3);
plotHandle4 = plot(sampData(5).massFracCorr238U235U, '.');
plotHandle4.MarkerSize = 14;
pH5 = subplot(5, 2, 4);
plotHandle5 = plot(sampData(7).massFracCorr238U235U, '.');
plotHandle5.MarkerSize = 14;
pH6 = subplot(5, 2, 5);
plotHandle6 = plot(sampData(9).massFracCorr238U235U, '.');
plotHandle6.MarkerSize = 14;
pH7 = subplot(5, 2, 6);
plotHandle7 = plot(sampData(11).massFracCorr238U235U, '.');
plotHandle7.MarkerSize = 14;
pH8 = subplot(5, 2, 7);
plotHandle8 = plot(sampData(13).massFracCorr238U235U, '.');
plotHandle8.MarkerSize = 14;
pH9 = subplot(5, 2, 8);
plotHandle9 = plot(sampData(15).massFracCorr238U235U, '.');
plotHandle9.MarkerSize = 14;
pH10 = subplot(5, 2, 9);
plotHandle10 = plot(sampData(17).massFracCorr238U235U, '.');
plotHandle10.MarkerSize = 14;
pH11 = subplot(5, 2, 10);
plotHandle11 = plot(sampData(33).massFracCorr238U235U, '.');
plotHandle11.MarkerSize = 14;


figHandle2 = figure('Name', '232Th229ThRatios');
pH12 = subplot(5, 2, 1);
plotHandle12 = plot(sampData(2).massFracCorr232Th229Th, '.');
plotHandle12.MarkerSize = 14;
plotHandle12.MarkerEdgeColor = 'r';
plotHandle12.MarkerFaceColor = 'r';
pH13 = subplot(5, 2, 2);
plotHandle13 = plot(sampData(4).massFracCorr232Th229Th, '.');
plotHandle13.MarkerSize = 14;
plotHandle13.MarkerEdgeColor = 'r';
plotHandle13.MarkerFaceColor = 'r';
pH14 = subplot(5, 2, 3);
plotHandle14 = plot(sampData(6).massFracCorr232Th229Th, '.');
plotHandle14.MarkerSize = 14;
plotHandle14.MarkerEdgeColor = 'r';
plotHandle14.MarkerFaceColor = 'r';
pH15 = subplot(5, 2, 4);
plotHandle15 = plot(sampData(8).massFracCorr232Th229Th, '.');
plotHandle15.MarkerSize = 14;
plotHandle15.MarkerEdgeColor = 'r';
plotHandle15.MarkerFaceColor = 'r';
pH16 = subplot(5, 2, 5);
plotHandle16 = plot(sampData(10).massFracCorr232Th229Th, '.');
plotHandle16.MarkerSize = 14;
plotHandle16.MarkerEdgeColor = 'r';
plotHandle16.MarkerFaceColor = 'r';
pH17 = subplot(5, 2, 6);
plotHandle17 = plot(sampData(12).massFracCorr232Th229Th, '.');
plotHandle17.MarkerSize = 14;
plotHandle17.MarkerEdgeColor = 'r';
plotHandle17.MarkerFaceColor = 'r';
pH18 = subplot(5, 2, 7);
plotHandle18 = plot(sampData(14).massFracCorr232Th229Th, '.');
plotHandle18.MarkerSize = 14;
plotHandle18.MarkerEdgeColor = 'r';
plotHandle18.MarkerFaceColor = 'r';
pH19 = subplot(5, 2, 8);
plotHandle19 = plot(sampData(16).massFracCorr232Th229Th, '.');
plotHandle19.MarkerSize = 14;
plotHandle19.MarkerEdgeColor = 'r';
plotHandle19.MarkerFaceColor = 'r';
pH20 = subplot(5, 2, 9);
plotHandle20 = plot(sampData(18).massFracCorr232Th229Th, '.');
plotHandle20.MarkerSize = 14;
plotHandle20.MarkerEdgeColor = 'r';
plotHandle20.MarkerFaceColor = 'r';
pH21 = subplot(5, 2, 10);
plotHandle21 = plot(sampData(34).massFracCorr232Th229Th, '.');
plotHandle21.MarkerSize = 14;
plotHandle21.MarkerEdgeColor = 'r';
plotHandle21.MarkerFaceColor = 'r';


figHandle3 = figure('Name', '229ThIntensities');
pH22 = subplot(5, 2, 1);
plotHandle22 = plot(sampData(2).OPZcorrInt(:, 1), '.');
plotHandle22.MarkerSize = 14;
plotHandle22.MarkerEdgeColor = 'm';
plotHandle22.MarkerFaceColor = 'm';
pH23 = subplot(5, 2, 2);
plotHandle23 = plot(sampData(4).OPZcorrInt(:, 1), '.');
plotHandle23.MarkerSize = 14;
plotHandle23.MarkerEdgeColor = 'm';
plotHandle23.MarkerFaceColor = 'm';
pH24 = subplot(5, 2, 3);
plotHandle24 = plot(sampData(6).OPZcorrInt(:, 1), '.');
plotHandle24.MarkerSize = 14;
plotHandle24.MarkerEdgeColor = 'm';
plotHandle24.MarkerFaceColor = 'm';
pH25 = subplot(5, 2, 4);
plotHandle25 = plot(sampData(8).OPZcorrInt(:, 1), '.');
plotHandle25.MarkerSize = 14;
plotHandle25.MarkerEdgeColor = 'm';
plotHandle25.MarkerFaceColor = 'm';
pH26 = subplot(5, 2, 5);
plotHandle26 = plot(sampData(10).OPZcorrInt(:, 1), '.');
plotHandle26.MarkerSize = 14;
plotHandle26.MarkerEdgeColor = 'm';
plotHandle26.MarkerFaceColor = 'm';
pH27 = subplot(5, 2, 6);
plotHandle27 = plot(sampData(12).OPZcorrInt(:, 1), '.');
plotHandle27.MarkerSize = 14;
plotHandle27.MarkerEdgeColor = 'm';
plotHandle27.MarkerFaceColor = 'm';
pH28 = subplot(5, 2, 7);
plotHandle28 = plot(sampData(14).OPZcorrInt(:, 1), '.');
plotHandle28.MarkerSize = 14;
plotHandle28.MarkerEdgeColor = 'm';
plotHandle28.MarkerFaceColor = 'm';
pH29 = subplot(5, 2, 8);
plotHandle29 = plot(sampData(16).OPZcorrInt(:, 1), '.');
plotHandle29.MarkerSize = 14;
plotHandle29.MarkerEdgeColor = 'm';
plotHandle29.MarkerFaceColor = 'm';
pH30 = subplot(5, 2, 9);
plotHandle30 = plot(sampData(18).OPZcorrInt(:, 1), '.');
plotHandle30.MarkerSize = 14;
plotHandle30.MarkerEdgeColor = 'm';
plotHandle30.MarkerFaceColor = 'm';
pH31 = subplot(5, 2, 10);
plotHandle31 = plot(sampData(34).OPZcorrInt(:, 1), '.');
plotHandle31.MarkerSize = 14;
plotHandle31.MarkerEdgeColor = 'm';
plotHandle31.MarkerFaceColor = 'm';


figHandle4 = figure('Name', '235UIntensities');
pH32 = subplot(5, 2, 1);
plotHandle32 = plot(sampData(1).OPZcorrInt(:, 3), '.');
plotHandle32.MarkerSize = 14;
plotHandle32.MarkerEdgeColor = 'k';
plotHandle32.MarkerFaceColor = 'k';
pH33 = subplot(5, 2, 2);
plotHandle33 = plot(sampData(3).OPZcorrInt(:, 3), '.');
plotHandle33.MarkerSize = 14;
plotHandle33.MarkerEdgeColor = 'k';
plotHandle33.MarkerFaceColor = 'k';
pH34 = subplot(5, 2, 3);
plotHandle34 = plot(sampData(5).OPZcorrInt(:, 3), '.');
plotHandle34.MarkerSize = 14;
plotHandle34.MarkerEdgeColor = 'k';
plotHandle34.MarkerFaceColor = 'k';
pH35 = subplot(5, 2, 4);
plotHandle35 = plot(sampData(7).OPZcorrInt(:, 3), '.');
plotHandle35.MarkerSize = 14;
plotHandle35.MarkerEdgeColor = 'k';
plotHandle35.MarkerFaceColor = 'k';
pH36 = subplot(5, 2, 5);
plotHandle36 = plot(sampData(9).OPZcorrInt(:, 3), '.');
plotHandle36.MarkerSize = 14;
plotHandle36.MarkerEdgeColor = 'k';
plotHandle36.MarkerFaceColor = 'k';
pH37 = subplot(5, 2, 6);
plotHandle37 = plot(sampData(11).OPZcorrInt(:, 3), '.');
plotHandle37.MarkerSize = 14;
plotHandle37.MarkerEdgeColor = 'k';
plotHandle37.MarkerFaceColor = 'k';
pH38 = subplot(5, 2, 7);
plotHandle38 = plot(sampData(13).OPZcorrInt(:, 3), '.');
plotHandle38.MarkerSize = 14;
plotHandle38.MarkerEdgeColor = 'k';
plotHandle38.MarkerFaceColor = 'k';
pH39 = subplot(5, 2, 8);
plotHandle39 = plot(sampData(15).OPZcorrInt(:, 3), '.');
plotHandle39.MarkerSize = 14;
plotHandle39.MarkerEdgeColor = 'k';
plotHandle39.MarkerFaceColor = 'k';
pH40 = subplot(5, 2, 9);
plotHandle40 = plot(sampData(17).OPZcorrInt(:, 3), '.');
plotHandle40.MarkerSize = 14;
plotHandle40.MarkerEdgeColor = 'k';
plotHandle40.MarkerFaceColor = 'k';
pH41 = subplot(5, 2, 10);
plotHandle41 = plot(sampData(33).OPZcorrInt(:, 3), '.');
plotHandle41.MarkerSize = 14;
plotHandle41.MarkerEdgeColor = 'k';
plotHandle41.MarkerFaceColor = 'k';

