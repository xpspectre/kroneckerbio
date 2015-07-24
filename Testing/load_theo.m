function data = load_theo
% Load and preprocess data (normalize from NONMEM form into timeseries/patient)
%% Load from csv
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

formatSpec = '%f%f%f%f%f%[^\n\r]';

fileID = fopen('theo.csv','r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

Theoph = table(dataArray{1:end-1}, 'VariableNames', {'Subject','Wt','Dose','Time','Conc'});

%% Convert to convenient form
ids = unique(Theoph.Subject);
nIds = length(ids);
dataOutCell = cell(nIds,4);
for i = 1:nIds
    id = ids(i);
    
    data = Theoph(Theoph.Subject == id,:);
    
    Wt = data.Wt(1);
    DosePer = data.Dose(1);
    Dose = Wt*DosePer;
    
    Time = data.Time;
    Conc = data.Conc;

    dataOutCell(i,:) = {id, Time, Conc, Dose};
end
data = cell2table(dataOutCell, 'VariableNames', {'ID','Time','Conc','Dose'});