%% Parse a folder of Element2 data files, perform on peak zeros subtraction
%  and calculate isotope ratios.
%  For 229Th-232Th-235U-238U runs

%% Setup for parsing text file
delimiter = '\t';
startRow = 6;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';


%% Sort files out of folder

datafolderstring = './data/112020_AC_MoreAcidBlanks/';
% don't forget to put a / at the end of the folder name!

fileStruct = dir([datafolderstring '*.txt']);
n.files = numel(fileStruct);
% minus one for the .txt file with the sample list generated from sequence
[temp.min, temp.indx] = min([fileStruct.bytes]); % sample list is smallest file
% fileStruct(temp.indx) = []; % delete sample list

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
isBlank = ~logical(1:length(data))';
isNotRinse = ~logical(1:length(data))';

for iRinse = 1:length(data)
    
    % Remove rinses
    temp.nameString1 = data(iRinse).fileName(4:end); 
    if strcmp(temp.nameString1, 'Rinse')
        isNotRinse(iRinse) = 0;
    end
    if strcmp(temp.nameString1, 'Wash')
        isNotRinse(iRinse) = 1;
    end
    temp.nameString1 = data(iRinse).fileName(7:end); 
    if strcmp(temp.nameString1, 'std')
        isNotRinse(iRinse) = 1;
    end
    
    temp.nameString1 = data(iRinse).fileName(7); 
    if strcmp(temp.nameString1, 'M')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, 'C')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '5')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '0')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '3')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, '2')
        isNotRinse(iRinse) = 1;
    end
    if strcmp(temp.nameString1, 'R')
        isNotRinse(iRinse) = 1;
    end
    
end

for iBlank = 1: length(data)
    
      % Flag washes to get backgrounds
    temp.nameString1 = data(iBlank).fileName(7); 
    if strcmp(temp.nameString1, 'M')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, 'C')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, '5')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, '0')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, '3')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, '2')
        isBlank(iBlank) = 1;
    end
    if strcmp(temp.nameString1, 'R')
        isBlank(iBlank) = 1;
    end
    
end


data = data(isNotRinse);
isWash = ~logical(1:length(data))';


for iWash = 1: length(data)
    
      % Flag washes to get backgrounds
    temp.sampNameString = data(iWash).fileName(4:end);
    if strcmp(temp.sampNameString, 'Wash')
        isWash(iWash) = 1;
    end
    
end

WashData = data(isWash);

for iBackgrounds = 1: length(WashData)
    
    % Find average background intensities
    BG232Th = WashData(iBackgrounds).intensities(:, 6);
    BG229Th = WashData(iBackgrounds).intensities(:, 9);
    BG238U = WashData(iBackgrounds).intensities(:, 8);
    BG235U = WashData(iBackgrounds).intensities(:, 7);
    
    WashData(iBackgrounds).AvgBG232Th = mean(BG232Th);
    WashData(iBackgrounds).AvgBG229Th = mean(BG229Th);
    WashData(iBackgrounds).AvgBG235U = mean(BG235U);
    WashData(iBackgrounds).AvgBG238U = mean(BG238U);
    
    WashData(iBackgrounds).BG232ThError = std(BG232Th)./ sqrt(length(BG232Th));
    WashData(iBackgrounds).BG229ThError = std(BG229Th)./ sqrt(length(BG229Th));
    WashData(iBackgrounds).BG235UError = std(BG235U)./ sqrt(length(BG235U));
    WashData(iBackgrounds).BG238UError = std(BG238U)./ sqrt(length(BG238U));
    
end


for iRun = 2:2:length(data)
    
    % OPZ subtraction
    data(iRun).OPZcorrInt = data(iRun).intensities(10:(end-10),:) - data(iRun-1).intensities(10:(end-10),:);
    
    % Flag standards
    temp.nameString = data(iRun).fileName(7:end); 
    if strcmp(temp.nameString, 'std')
        isSTD(iRun) = 1;
    end
    
end
    
BlankData = data(isBlank);

%% calculate ratios and average ratios


for iRatio = 2:2:length(data)

    Th232 = (data(iRatio).OPZcorrInt(:,6));
    Th229 = (data(iRatio).OPZcorrInt(:,9));
    data(iRatio).Th232Th229Ratio = Th232 ./ Th229 ;
    
    
    U235 = (data(iRatio).OPZcorrInt(:,7));
    U238 = (data(iRatio).OPZcorrInt(:,8));
    data(iRatio).U238U235Ratio = U238 ./ U235 ;
    
    data(iRatio).U235MeanIntensity = mean(U235);
    data(iRatio).Th229MeanIntensity = mean(Th229);
    data(iRatio).U238MeanIntensity = mean(U238);
    data(iRatio).Th232MeanIntensity = mean(Th232);
    
    data(iRatio).avg232Th229ThRatios = mean(data(iRatio).Th232Th229Ratio);
    data(iRatio).avg238U235URatios = mean(data(iRatio).U238U235Ratio);
end

%% Calculate mean blank intensities

for iBlanks = 1: length(BlankData)
    
    BlankData(iBlanks).OPZcorrInt(BlankData(iBlanks).OPZcorrInt < 0)  = 0;
    
    BlankTh229 = BlankData(iBlanks).intensities(:, 9);
    BlankTh232 = BlankData(iBlanks).intensities(:, 6);
    BlankU235 = BlankData(iBlanks).intensities(:, 7);
    BlankU238 = BlankData(iBlanks).intensities(:, 8);
    
    BlankData(iBlanks).U235MeanIntensity = mean(BlankU235);
    BlankData(iBlanks).Th229MeanIntensity = mean(BlankTh229);
    BlankData(iBlanks).U238MeanIntensity = mean(BlankU238);
    BlankData(iBlanks).Th232MeanIntensity = mean(BlankTh232);
    
end

%% Calculate betas and correct average ratios for mass fractionation

CRM112aURatio = 137.844;
SRMU970URatio = 186.78;
mass238 = 238;
mass235 = 235;
stdData = data(isSTD);


for iStd = 1:length(stdData)
    
    stdData(iStd).Th232U238 = stdData(iStd).OPZcorrInt(:, 6) ./ stdData(iStd).OPZcorrInt(:, 8);
    stdData(iStd).avgTh232U238 = mean(stdData(iStd).Th232U238);
    stdData(iStd).U238U235StdError = std(stdData(iStd).U238U235Ratio) / (sqrt(length(stdData(iStd).U238U235Ratio)));
    stdData(iStd).U238U235FracUncert = mean(stdData(iStd).U238U235StdError ./ stdData(iStd).U238U235Ratio);
    stdData(iStd).Beta = (log(CRM112aURatio) - log(stdData(iStd).U238U235Ratio)) ./ log(mass235/mass238);
    stdData(iStd).BetaAvg = mean(stdData(iStd).Beta);
    stdData(iStd).BetaStdError = std(stdData(iStd).Beta) / sqrt(length(stdData(iStd).Beta));

end

%% Calculate U and Th blanks

StdUconc = 1.224 / 5;
StdThconc = 2.868 / 5;



%% Plotting backgrounds to check for outliers
% 
% figHandle = figure('Name', 'Backgrounds');
% pH2 = subplot(3, 2, 1);
% hold on
% plotHandle2 = plot(WashData(1).intensities(:,9), 'color', '#0072BD');
% plotHandle2 = plot(WashData(1).intensities(:,6), 'color', '#EDB120');
% plotHandle2 = plot(WashData(1).intensities(:,7), 'color', '#D95319');
% plotHandle2 = plot(WashData(1).intensities(:,8), '.', 'color', 'black');
% title('Wash-01');
% plotHandle2.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH3 = subplot(3, 2, 2);
% hold on
% plotHandle3 = plot(WashData(3).intensities(:,9), 'color', '#0072BD');
% plotHandle3 = plot(WashData(3).intensities(:,6), 'color', '#EDB120');
% plotHandle3 = plot(WashData(3).intensities(:,7), 'color', '#D95319');
% plotHandle3 = plot(WashData(3).intensities(:,8), '.', 'color', 'black');
% title('Wash-05');
% plotHandle3.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH4 = subplot(3, 2, 3);
% hold on
% plotHandle4 = plot(WashData(9).intensities(:,9), 'color', '#0072BD');
% plotHandle4 = plot(WashData(9).intensities(:,6), 'color', '#EDB120');
% plotHandle4 = plot(WashData(9).intensities(:,7), 'color', '#D95319');
% plotHandle4 = plot(WashData(9).intensities(:,8), '.', 'color', 'black');
% title('Wash-17');
% plotHandle4.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH5 = subplot(3, 2, 4);
% hold on
% plotHandle5 = plot(WashData(15).intensities(:,9), 'color', '#0072BD');
% plotHandle5 = plot(WashData(15).intensities(:,6), 'color', '#EDB120');
% plotHandle5 = plot(WashData(15).intensities(:,7), 'color', '#D95319');
% plotHandle5 = plot(WashData(15).intensities(:,8), '.', 'color', 'black');
% title('Wash-29');
% plotHandle5.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH6 = subplot(3, 2, 5);
% hold on
% plotHandle6 = plot(WashData(21).intensities(:,9), 'color', '#0072BD');
% plotHandle6 = plot(WashData(21).intensities(:,6), 'color', '#EDB120');
% plotHandle6 = plot(WashData(21).intensities(:,7), 'color', '#D95319');
% plotHandle6 = plot(WashData(21).intensities(:,8), '.', 'color', 'black');
% title('Wash-41');
% plotHandle6.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH7 = subplot(3, 2, 6);
% hold on
% plotHandle7 = plot(WashData(26).intensities(:,9), 'color', '#0072BD');
% plotHandle7 = plot(WashData(26).intensities(:,6), 'color', '#EDB120');
% plotHandle7 = plot(WashData(26).intensities(:,7), 'color', '#D95319');
% plotHandle7 = plot(WashData(26).intensities(:,8), '.', 'color', 'black');
% title('Wash-51');
% plotHandle7.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% 
% 
% %% Plotting blank data to search for outliers
% 
% figHandle2 = figure('Name', 'Blanks');
% pH8 = subplot(3, 2, 1);
% hold on
% plotHandle2 = plot(WashData(1).intensities(:,9), 'color', '#0072BD');
% plotHandle2 = plot(WashData(1).intensities(:,6), 'color', '#EDB120');
% plotHandle2 = plot(WashData(1).intensities(:,7), 'color', '#D95319');
% plotHandle2 = plot(WashData(1).intensities(:,8), '.', 'color', 'black');
% title('Wash-01');
% plotHandle2.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH9 = subplot(3, 2, 2);
% hold on
% plotHandle3 = plot(WashData(3).intensities(:,9), 'color', '#0072BD');
% plotHandle3 = plot(WashData(3).intensities(:,6), 'color', '#EDB120');
% plotHandle3 = plot(WashData(3).intensities(:,7), 'color', '#D95319');
% plotHandle3 = plot(WashData(3).intensities(:,8), '.', 'color', 'black');
% title('Wash-05');
% plotHandle3.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH4 = subplot(3, 2, 3);
% hold on
% plotHandle4 = plot(WashData(9).intensities(:,9), 'color', '#0072BD');
% plotHandle4 = plot(WashData(9).intensities(:,6), 'color', '#EDB120');
% plotHandle4 = plot(WashData(9).intensities(:,7), 'color', '#D95319');
% plotHandle4 = plot(WashData(9).intensities(:,8), '.', 'color', 'black');
% title('Wash-17');
% plotHandle4.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH5 = subplot(3, 2, 4);
% hold on
% plotHandle5 = plot(WashData(15).intensities(:,9), 'color', '#0072BD');
% plotHandle5 = plot(WashData(15).intensities(:,6), 'color', '#EDB120');
% plotHandle5 = plot(WashData(15).intensities(:,7), 'color', '#D95319');
% plotHandle5 = plot(WashData(15).intensities(:,8), '.', 'color', 'black');
% title('Wash-29');
% plotHandle5.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH6 = subplot(3, 2, 5);
% hold on
% plotHandle6 = plot(WashData(21).intensities(:,9), 'color', '#0072BD');
% plotHandle6 = plot(WashData(21).intensities(:,6), 'color', '#EDB120');
% plotHandle6 = plot(WashData(21).intensities(:,7), 'color', '#D95319');
% plotHandle6 = plot(WashData(21).intensities(:,8), '.', 'color', 'black');
% title('Wash-41');
% plotHandle6.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% pH7 = subplot(3, 2, 6);
% hold on
% plotHandle7 = plot(WashData(26).intensities(:,9), 'color', '#0072BD');
% plotHandle7 = plot(WashData(26).intensities(:,6), 'color', '#EDB120');
% plotHandle7 = plot(WashData(26).intensities(:,7), 'color', '#D95319');
% plotHandle7 = plot(WashData(26).intensities(:,8), '.', 'color', 'black');
% title('Wash-51');
% plotHandle7.LineWidth = 1;
% ylim ([0 10000]);
% hold off
% 
