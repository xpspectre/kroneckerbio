function str = cellstr2str(cellstr)
% Convert cell array of strings to single string with cells as 'cell1, cell2, etc.'

cellstr = vec(cellstr);
str = cellfun(@(x) [x ', '], cellstr, 'UniformOutput', false);
str = cell2mat(str');
str(end-1:end) = []; % trim off trailing ', '